"""
    read_nc_variable(dataset, path::String; bounds=false)

Read a variable from a NetCDF dataset, handling nested groups via `/` separators.

If `bounds=true`, reshape the result so that footprint vertices are in the last dimension,
returning an N×4 matrix (N soundings, 4 corners).
"""
function read_nc_variable(dataset, path::String; bounds::Bool=false)
    var = _nc_variable(dataset, path)
    bounds && return _read_bounds(var)
    data = var.var[:]
    occursin('/', path) ? vec(data) : data
end

# Resolve paths once for raw reads, CF reads, and variable attributes.
function _nc_variable(dataset, path::String)
    parts = split(strip(path, '/'), '/')
    node = dataset
    for group in parts[1:end-1]
        node = node.group[group]
    end
    node[parts[end]]
end

function _read_bounds(var)
    shape = size(var)
    axes4 = findall(==(4), collect(shape))
    isempty(axes4) && return var.var[:]
    # Size alone is ambiguous for four soundings. Prefer a named vertex axis;
    # retain the historical first size-4 axis when dimension names give no clue.
    names = lowercase.(String.(dimnames(var)))
    vertex_axes = filter(axes4) do axis
        name = names[axis]
        occursin(r"vertex|vertices|corner|bound|bnd|^n?v$|^nv4$", name)
    end
    axis = length(vertex_axes) == 1 ? only(vertex_axes) : first(axes4)
    data = reshape(var.var[:], shape)
    axis == ndims(var) && return reshape(data, :, 4)
    axis == 1 && return transpose(reshape(data, 4, :))
    order = (filter(!=(axis), collect(1:ndims(var)))..., axis)
    reshape(permutedims(data, order), :, 4)
end

_read_flat_variable(dataset, varname::String, bounds::Bool) =
    bounds ? _read_bounds(dataset[varname]) : dataset[varname].var[:]

"""
    read_variable_from_file(path, variable; bounds=false)

Read a variable from a NetCDF-like file. Falls back to the system `libnetcdf`
for flat HDF4/HDF-EOS variables when the Julia NetCDF build cannot open them.
"""
function read_variable_from_file(path::String, variable::String; bounds::Bool=false)
    if _looks_like_hdf4(path) && !bounds
        return read_system_netcdf_variable(path, variable)
    end
    try
        ds = Dataset(path)
        try
            return read_nc_variable(ds, variable; bounds=bounds)
        finally
            close(ds)
        end
    catch e
        if bounds
            rethrow(e)
        end
        try
            return read_system_netcdf_variable(path, variable)
        catch fallback_error
            error("Could not read variable '$variable' from '$path' with NCDatasets ($(sprint(showerror, e))) or system libnetcdf fallback ($(sprint(showerror, fallback_error)))")
        end
    end
end

# HDF4 (.hdf) files are not supported by the bundled NetCDF_jll on most
# platforms, so an `NCDatasets.Dataset(path)` attempt fails after a 150–500 ms
# NFS round-trip. Skip that try and go straight to the system libnetcdf path,
# which links against an HDF4-capable build.
_looks_like_hdf4(path::AbstractString) = endswith(lowercase(path), ".hdf")

"""
    read_array_from_file(path, variable)

Read a variable while preserving its native dimensionality.
"""
function read_array_from_file(path::String, variable::String)
    if _looks_like_hdf4(path)
        return read_system_netcdf_variable(path, variable)
    end
    try
        ds = Dataset(path)
        try
            return _read_nc_array(ds, variable)
        finally
            close(ds)
        end
    catch e
        try
            return read_system_netcdf_variable(path, variable)
        catch fallback_error
            error("Could not read variable '$variable' from '$path' with NCDatasets ($(sprint(showerror, e))) or system libnetcdf fallback ($(sprint(showerror, fallback_error)))")
        end
    end
end

"""
    read_arrays_from_file(path, variables) -> Vector

Read multiple variables from a single file open. Each variable is read in full
(same semantics as `read_array_from_file`). Falls back to per-variable reads
through the system libnetcdf path if NCDatasets cannot open the file.
"""
function read_arrays_from_file(path::String,
                               variables::AbstractVector{<:AbstractString})
    if _looks_like_hdf4(path)
        return read_system_netcdf_variables(path, variables)
    end
    try
        ds = Dataset(path)
        try
            return Any[_read_nc_array(ds, String(v)) for v in variables]
        finally
            close(ds)
        end
    catch e
        try
            return read_system_netcdf_variables(path, variables)
        catch fallback_error
            error("Could not read variables from '$path' with NCDatasets ($(sprint(showerror, e))) or system libnetcdf fallback ($(sprint(showerror, fallback_error)))")
        end
    end
end

_read_nc_array(dataset, path::String) = _read_full_array(_nc_variable(dataset, path))

function _read_full_array(var)
    selectors = ntuple(_ -> Colon(), ndims(var))
    var[selectors...]
end

# Science inputs must honor CF packing and missing-value attributes. The public
# read_nc_variable helper retains its raw-read semantics for existing callers.
function _read_l2_values(dataset, path::String, ::Type{T}) where {T}
    data = _read_nc_array(dataset, path)
    map(x -> ismissing(x) ? T(NaN) : T(x), vec(data))
end

"""
    read_nc_attribute(dataset, path::String, attr::String)

Read an attribute from a possibly group-nested variable.
"""
read_nc_attribute(dataset, path::String, attr::String) =
    _nc_variable(dataset, path).attrib[attr]


# Define the per-snapshot averaging-window bounds along the unlimited `time` dim.
# `time` itself holds the snapshot reference (start of interval, or window center for a
# rolling mean); `window_start`/`window_end` make the actual averaging interval explicit
# and self-documenting. Both are `(time,)` so they concatenate with ncrcat like `time`.
function _define_window_bounds!(ds, n_times::Int)
    ws = defVar(ds, "window_start", Float32, ("time",),
                attrib=["units" => "days since 1970-01-01",
                         "long_name" => "Averaging window start (inclusive)"])
    we = defVar(ds, "window_end", Float32, ("time",),
                attrib=["units" => "days since 1970-01-01",
                         "long_name" => "Averaging window end (inclusive)"])
    ws[1:n_times] = zeros(Float32, n_times)
    we[1:n_times] = zeros(Float32, n_times)
    nothing
end

# Transfer ownership to the caller only after all definitions succeed.
function _with_output_setup(f, outfile::String)
    ds = Dataset(outfile, "c")
    try
        f(ds)
    catch
        close(ds)
        rethrow()
    end
end

function _close_output_dataset!(ds, metadata::InputMetadata)
    ds === nothing && return nothing
    try
        _write_input_metadata!(ds, metadata)
    finally
        close(ds)
    end
    nothing
end

"""
    create_output_dataset(outfile, grid_spec, n_times, grid_vars, compute_std) -> Dataset, Dict

Create the output NetCDF4 file with dimensions, coordinate variables, and data variable
definitions. Returns the dataset and a Dict mapping variable names to NCDatasets variables.
The global `command` and `filters` attributes record the invocation and configured
quality-filter rules when provided.
"""
function create_output_dataset(outfile::String, grid_spec::RectangularGridSpec,
                               n_times::Int, grid_vars::OrderedDict{String,String},
                               compute_std::Bool;
                               window_bounds::Bool=false,
                               time_long_name::AbstractString="Time (UTC), start of interval",
                               command::Union{Nothing,String}=nothing,
                               filters::Vector{FilterRule}=FilterRule[],
                               output_options::NetCDFOptions=NetCDFOptions())
    _create_output_dataset(outfile, grid_spec, n_times, grid_vars, compute_std;
        window_bounds=window_bounds, time_long_name=time_long_name,
        command=command, filters=filters, output_options=output_options)
end

function create_output_dataset(outfile::String, grid_spec::CubedSphereGridSpec,
                               n_times::Int, grid_vars::OrderedDict{String,String},
                               compute_std::Bool;
                               window_bounds::Bool=false,
                               time_long_name::AbstractString="Time (UTC), start of interval",
                               command::Union{Nothing,String}=nothing,
                               filters::Vector{FilterRule}=FilterRule[],
                               output_options::NetCDFOptions=NetCDFOptions())
    _create_output_dataset(outfile, grid_spec, n_times, grid_vars, compute_std;
        window_bounds=window_bounds, time_long_name=time_long_name,
        command=command, filters=filters, output_options=output_options)
end

_output_title(::RectangularGridSpec) = "Gridded satellite data"
_output_title(::CubedSphereGridSpec) = "Gridded satellite data (cubed sphere)"
_payload_attributes(::RectangularGridSpec) = Pair{String,String}[]
_payload_attributes(::CubedSphereGridSpec) =
    ["coordinates" => "lons lats", "grid_mapping" => "cubed_sphere"]

function _create_output_dataset(outfile, spec, n_times, grid_vars, compute_std;
                                 window_bounds, time_long_name, command, filters, output_options)
    storage = _storage_kwargs(output_options, spec)
    _with_output_setup(outfile) do ds
        write_coordinates!(ds, spec)
        # Unlimited time supports concatenation with ncrcat. Placeholder times
        # size all dependent variables; the gridders write the actual dates.
        defDim(ds, "time", Inf)
        time = defVar(ds, "time", Float32, ("time",);
                      attrib=["units" => "days since 1970-01-01", "long_name" => time_long_name])
        time[1:n_times] = zeros(Float32, n_times)
        window_bounds && _define_window_bounds!(ds, n_times)
        ds.attrib["title"] = _output_title(spec)
        ds.attrib["created_with"] = "SatelliteGridding.jl"
        _set_provenance_attributes!(ds, command, filters)

        dims = (payload_dims(spec)..., "time")
        attributes = _payload_attributes(spec)
        nc_vars = Dict{String,Any}()
        for key in keys(grid_vars)
            nc_vars[key] = defVar(ds, key, Float32, dims; storage..., attrib=attributes)
            if compute_std
                nc_vars[key * "_std"] = defVar(ds, key * "_std", Float32, dims;
                                               storage..., attrib=attributes)
            end
        end
        nc_vars["n"] = defVar(ds, "n", Float32, dims; storage...,
            attrib=vcat(["units" => "", "long_name" => "Number of pixels in average"], attributes))
        ds, nc_vars
    end
end

function _set_provenance_attributes!(ds, command::Union{Nothing,String},
                                     filters::Vector{FilterRule})
    for (key, value) in _software_provenance()
        ds.attrib[key] = value
    end
    command === nothing || (ds.attrib["command"] = command)
    ds.attrib["filters"] = isempty(filters) ? "none" :
        join((_format_filter_rule(rule) for rule in filters), "; ")
    nothing
end

function _format_filter_rule(rule::FilterRule)
    if rule.op === :between
        "$(rule.variable): $(rule.lo) < x < $(rule.hi)"
    else
        operator = rule.op === :lt ? "<" : rule.op === :gt ? ">" : "=="
        "$(rule.variable): x $operator $(rule.lo)"
    end
end

"""
    copy_variable_attributes!(nc_vars, dataset, grid_vars)

Copy attributes (units, long_name, etc.) from the first input file to the output variables.
"""
function copy_variable_attributes!(nc_vars::Dict, dataset, grid_vars::OrderedDict{String,String};
                                   attributes=DEFAULT_VARIABLE_ATTRIBUTES)
    for (key, value) in grid_vars
        for at in attributes
            try
                nc_vars[key].attrib[at] = read_nc_attribute(dataset, value, at)
            catch
            end
        end
    end
end
