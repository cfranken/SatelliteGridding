# Only small metadata variables are read; never load per-sounding version arrays.
const DEFAULT_VERSION_SOURCES = [
    "@product_version", "@ProductVersion", "@version", "@VersionID",
    "@build_id", "@BuildId", "@BuildID", "@CollectionVersion",
    "Metadata/BuildId", "@gesdisc_collection", "Metadata/CollectionLabel", "@CollectionLabel",
]
const DEFAULT_VARIABLE_ATTRIBUTES = [
    "units", "long_name", "valid_range", "description", "unit", "longname",
]

function _validate_metadata_options(settings)
    settings isa AbstractDict || throw(ArgumentError("metadata must be a table"))
    for (key, value) in settings
        key in ("version_sources", "attributes", "variable_attributes") ||
            throw(ArgumentError("Unknown [metadata] option: $key"))
        value isa AbstractVector && all(x -> x isa AbstractString, value) ||
            throw(ArgumentError("metadata.$key must be an array of strings"))
    end
    # These describe the source's storage, not the gridded Float32 values.
    for name in get(settings, "variable_attributes", DEFAULT_VARIABLE_ATTRIBUTES)
        name in ("_FillValue", "missing_value", "scale_factor", "add_offset", "_Unsigned") &&
            throw(ArgumentError("Cannot copy storage attribute '$name' into gridded variables"))
    end
    nothing
end

mutable struct InputMetadata
    settings::Dict{String,Any}
    versions::Dict{String,Vector{Any}}
    attributes::Dict{String,Vector{Any}}
    product_versions::Set{String}
    warned::Set{String}
end

function InputMetadata(config::DataSourceConfig)
    settings = Dict{String,Any}(get(config.options, "metadata", Dict()))
    _validate_metadata_options(settings)
    InputMetadata(settings, Dict{String,Vector{Any}}(), Dict{String,Vector{Any}}(),
                  Set{String}(), Set{String}())
end

function _git_output(root, args...)
    try
        # Avoid accidentally describing an enclosing repository in source archives.
        top = strip(read(pipeline(`git -C $root rev-parse --show-toplevel`; stderr=devnull), String))
        realpath(top) == realpath(root) || return "unknown"
        strip(read(pipeline(`git -C $root $args`; stderr=devnull), String))
    catch
        "unknown"
    end
end

function _software_provenance(; root=dirname(@__DIR__))
    project = TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
    tag = _git_output(root, "describe", "--tags", "--abbrev=0")
    commit = _git_output(root, "rev-parse", "HEAD")
    state = _git_output(root, "status", "--porcelain", "--untracked-files=no")
    Dict("software_version" => string(project["version"]),
         "git_tag" => tag, "git_commit" => commit,
         "git_dirty" => state == "unknown" ? "unknown" : string(!isempty(state)),
         "julia_version" => string(VERSION))
end

# Selectors: @attribute, Group/@attribute, or Group/small_variable.
function _metadata_value(ds, selector::AbstractString)
    parts = split(strip(selector, '/'), '/')
    node = ds
    for group in parts[1:end-1]
        haskey(node.group, group) || return nothing
        node = node.group[group]
    end
    name = parts[end]
    if startswith(name, "@")
        attr = name[2:end]
        return haskey(node.attrib, attr) ? node.attrib[attr] : nothing
    end
    haskey(node, name) || return nothing
    var = node[name]
    length(var) <= 1024 || throw(ArgumentError("Metadata variable '$selector' has more than 1024 elements"))
    data = ndims(var) == 0 ? var[] : var[ntuple(_ -> Colon(), ndims(var))...]
    data isa AbstractArray{Char} && return strip(join(vec(data)), ['\0', ' '])
    data isa AbstractArray && length(data) == 1 && return only(data)
    data
end

function _capture_input_metadata!(metadata::InputMetadata, ds)
    found_version = false
    for (dest, selectors) in (
        (metadata.versions, get(metadata.settings, "version_sources", DEFAULT_VERSION_SOURCES)),
        (metadata.attributes, get(metadata.settings, "attributes", String[])),
    )
        for selector in selectors
            try
                value = _metadata_value(ds, selector)
                (value === nothing || ismissing(value)) && continue
                value isa AbstractString && (value = strip(value, ['\0', ' ']))
                value == "" && continue
                values = get!(dest, selector, Any[])
                any(v -> isequal(v, value), values) || push!(values, value)
                if dest === metadata.versions && !found_version
                    push!(metadata.product_versions, string(value))
                    found_version = true
                end
            catch e
                if !(selector in metadata.warned)
                    @warn "Could not read optional input metadata" selector error=sprint(showerror, e)
                    push!(metadata.warned, selector)
                end
            end
        end
    end
    nothing
end
_capture_input_metadata!(::Nothing, ds) = nothing

function _write_input_metadata!(ds, metadata::InputMetadata)
    for (key, value) in _input_metadata_attributes(metadata)
        ds.attrib[key] = value
    end
    nothing
end

function _input_metadata_attributes(metadata::InputMetadata)
    versions = sort!(collect(metadata.product_versions))
    Dict("source_versions" => JSON.json(metadata.versions),
         "source_metadata" => JSON.json(metadata.attributes),
         "satellite_product_version" => isempty(versions) ? "unknown" : join(versions, "; "))
end
