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

function _read_bounds(var, raw=var.var[:])
    shape = size(var)
    axes4 = findall(==(4), collect(shape))
    isempty(axes4) && return raw
    # Size alone is ambiguous for four soundings. Prefer a named vertex axis;
    # retain the historical first size-4 axis when dimension names give no clue.
    names = lowercase.(String.(dimnames(var)))
    vertex_axes = filter(axes4) do axis
        name = names[axis]
        occursin(r"vertex|vertices|corner|bound|bnd|^n?v$|^nv4$", name)
    end
    axis = length(vertex_axes) == 1 ? only(vertex_axes) : first(axes4)
    data = reshape(raw, shape)
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

# Preserve the original coordinate values and precision for unpacked data.
function _mask_geolocation_missing(raw::Vector{T}, is_missing) where {T<:AbstractFloat}
    replace!(x -> is_missing(x) ? T(NaN) : x, raw)
end

function _mask_geolocation_missing(raw::Vector{T}, is_missing) where {T<:Integer}
    any(is_missing, raw) || return raw
    # Keep valid unpacked integers as integers so they participate in arithmetic
    # with the same precision as in a file without any missing coordinates.
    masked = Vector{Union{T,Float64}}(undef, length(raw))
    @inbounds for i in eachindex(raw)
        masked[i] = is_missing(raw[i]) ? NaN : raw[i]
    end
    return masked
end

_decode_geolocation(raw::Vector, is_missing, ::Nothing, ::Nothing) =
    _mask_geolocation_missing(raw, is_missing)

function _decode_geolocation(raw::Vector{T}, is_missing, scale, offset) where {T}
    scale_type = isnothing(scale) ? T : typeof(scale)
    offset_type = isnothing(offset) ? T : typeof(offset)
    decoded_type = float(promote_type(T, scale_type, offset_type))
    return map(raw) do x
        is_missing(x) && return decoded_type(NaN)
        scaled = isnothing(scale) ? x : x * scale
        decoded_type(isnothing(offset) ? scaled : scaled + offset)
    end
end

"""
    _read_l2_geolocation(dataset, path; bounds=false)

Read L2 coordinates in physical degrees, masking storage fill and missing
codes before applying any CF scale factor or offset. With no declared
`_FillValue`, use the NetCDF default for the storage type. Unpacked floating
point coordinates keep their original precision.
"""
function _read_l2_geolocation(dataset, path::String; bounds::Bool=false)
    var = _nc_variable(dataset, path)
    stored = var.var[:]
    raw = stored isa Number ? [stored] : vec(stored)
    fill = get(var.attrib, "_FillValue", NCDatasets.fillvalue(eltype(raw)))
    missing_values = get(var.attrib, "missing_value", ())
    scale = get(var.attrib, "scale_factor", nothing)
    offset = get(var.attrib, "add_offset", nothing)
    # NaN already represents invalid coordinates. Numeric equality also
    # masks both signed-zero representations of a zero fill code.
    is_missing(x) = x == fill || any(==(x), missing_values)
    data = _decode_geolocation(raw, is_missing, scale, offset)
    return bounds ? _read_bounds(var, data) : data
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
function create_output_dataset(outfile::String, grid_spec::AbstractGridSpec,
                               n_times::Int, grid_vars::OrderedDict{String,String},
                               compute_std::Bool;
                               window_bounds::Bool=false,
                               time_long_name::AbstractString="Time (UTC), start of interval",
                               command::Union{Nothing,String}=nothing,
                               filters::Vector{FilterRule}=FilterRule[],
                               output_options::NetCDFOptions=NetCDFOptions(),
                               spherical_ka_std::Bool=false)
    _create_output_dataset(outfile, grid_spec, n_times, grid_vars, compute_std;
        window_bounds=window_bounds, time_long_name=time_long_name,
        command=command, filters=filters, output_options=output_options, spherical_ka_std)
end

_output_title(::AbstractGridSpec) = "Gridded satellite data"
_output_title(::CubedSphereGridSpec) = "Gridded satellite data (cubed sphere)"
_output_title(::HEALPixGridSpec) = "Gridded satellite data (HEALPix)"
_payload_attributes(::AbstractGridSpec) = Pair{String,String}[]
_payload_attributes(::CubedSphereGridSpec) =
    ["coordinates" => "lons lats", "grid_mapping" => "cubed_sphere"]
_payload_attributes(::HEALPixGridSpec) =
    ["coordinates" => "lat lon cell_ids", "grid_mapping" => "healpix"]

_weight_long_name(::AbstractGridSpec) = "Number of pixels in average"
_weight_long_name(::HEALPixGridSpec) = "Sum of fractional footprint overlap weights"

function _create_output_dataset(outfile, spec, n_times, grid_vars, compute_std;
                                 window_bounds, time_long_name, command, filters, output_options,
                                 spherical_ka_std)
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
        rolling_comment = spherical_ka_std ?
            "Rolling windows combine centered daily moments. A non-finite variance or variance input in accumulator precision is written as fill." :
            "Rolling windows reconstruct variance from raw moments, where Float32 cancellation can produce a small nonzero value."
        single_comment = spherical_ka_std ?
            "Zero for a cell touched by a single footprint in tiled mode when its moments are finite. " :
            "Zero for a cell touched by a single footprint in tiled mode. "
        std_attributes = vcat(attributes, [
            "long_name" => "Weighted population standard deviation",
            "comment" => "sqrt(M2/n), where M2 is the weighted sum of squared deviations and n is the sum of weights. Computed without Bessel correction. " * single_comment * rolling_comment,
        ])
        nc_vars = Dict{String,Any}()
        for key in keys(grid_vars)
            nc_vars[key] = defVar(ds, key, Float32, dims; storage..., attrib=attributes)
            if compute_std
                nc_vars[key * "_std"] = defVar(ds, key * "_std", Float32, dims;
                                               storage..., attrib=std_attributes)
            end
        end
        nc_vars["n"] = defVar(ds, "n", Float32, dims; storage...,
            attrib=vcat(["units" => "1", "long_name" => _weight_long_name(spec)], attributes))
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

const _SOURCE_VARIABLE_REFERENCE_ATTRIBUTES = (
    "coordinates", "grid_mapping", "cell_measures", "bounds", "ancillary_variables", "formula_terms",
    "geometry", "node_coordinates", "climatology",
)

"""
    copy_variable_attributes!(nc_vars, dataset, grid_vars)

Copy attributes (units, long_name, etc.) from the first input file to the output variables.
Source variable references are skipped. If a standard-deviation variable exists,
it inherits the source units independently of the configured attribute list.
"""
function copy_variable_attributes!(nc_vars::Dict, dataset, grid_vars::OrderedDict{String,String};
                                   attributes=DEFAULT_VARIABLE_ATTRIBUTES)
    for (key, value) in grid_vars
        for at in attributes
            at in _SOURCE_VARIABLE_REFERENCE_ATTRIBUTES && continue
            try
                nc_vars[key].attrib[at] = read_nc_attribute(dataset, value, at)
            catch
            end
        end
        std_key = key * "_std"
        if haskey(nc_vars, std_key) && !haskey(grid_vars, std_key)
            try
                nc_vars[std_key].attrib["units"] = read_nc_attribute(dataset, value, "units")
            catch
            end
        end
    end
end
