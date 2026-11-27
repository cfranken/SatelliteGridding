"""
    grid_center(config, grid_spec, time_spec; kwargs...)

Grid data using center coordinates only (no footprint bounds), suitable for
MODIS-style data where each pixel maps to exactly one rectangular grid cell.

To archive Float64 center results exactly, construct `GridSpec` with Float64
bounds and spacings for Float64 accumulation and set `[output]`
`precision = "float64"` and `significant_digits_center = "none"`.
The center CLI uses a Float32 grid. A cell below `min_count` can have `n > 0`
alongside missing science values.

# Keyword Arguments
- `geo_table::Union{Nothing,String}=nothing`: Path to a legacy monolithic geolocation lookup table (NetCDF)
- `geo_cache::Union{Nothing,String}=nothing`: Directory for generated per-tile MODIS geolocation cache
- `geo_provider::Union{Symbol,String}=:auto`: `:auto`, `:variables`, `:lut`, or `:modis`
- `veg_indices::Bool=false`: Compute vegetation indices (EVI, NDVI, NIRv, NDWI)
- `verbose::Bool=true`: Log run status with `@info`.
- `outfile::String="gridded_output.nc"`: Output file path
- `keep_going::Bool=false`: Continue after per-file errors. By default, the
  first file error stops the run.
- `progress::Bool=stderr isa Base.TTY`: Draw progress bars (on stderr). By default they're
  drawn only when stderr is a terminal: redirected to a file or a pipe, as in batch jobs, every
  redraw would become new lines.
"""
function grid_center(config::DataSourceConfig, grid_spec::AbstractGridSpec{T}, time_spec::TimeSpec;
                     geo_table::Union{Nothing,String}=nothing,
                     geo_cache::Union{Nothing,String}=nothing,
                     geo_provider::Union{Symbol,String}=:auto,
                     veg_indices::Bool=false,
                     compute_std::Bool=false,
                     outfile::String="gridded_output.nc",
                     keep_going::Bool=false,
                     progress::Bool=_default_progress(),
                     verbose::Bool=true,
                     command::Union{Nothing,String}=nothing) where {T}
    _grid_center(config, grid_spec, time_spec; geo_table, geo_cache, geo_provider, veg_indices,
                 compute_std, outfile, keep_going, progress, verbose, command)
end

function _grid_center(config, spec, time_spec; kwargs...)
    error("Center-coordinate gridding requires a rectangular lat/lon grid. Got $(typeof(spec)).")
end

const _VEGETATION_NAMES = ("EVI", "NDVI", "NIRv", "NDWI")

"""Coordinates and mapped indices for the last tile returned by a center provider."""
struct _CenterIndexCache{I}
    # Used only for identity checks. A provider can return CF-decoded arrays
    # on a tile miss and its stored, converted arrays on a cache hit.
    lat::AbstractArray
    lon::AbstractArray
    idx::I
    lat_idx::Vector{Int32}
    lon_idx::Vector{Int32}
end

_center_cache(provider) = Ref{Union{Nothing,_CenterIndexCache}}(nothing)
_center_cache(::ModisSinusoidalGeolocation) =
    Ref{Union{Nothing,_CenterIndexCache{Vector{CartesianIndex{2}}}}}(nothing)

function _grid_center(config, grid_spec::RectangularGridSpec{T}, time_spec;
                      geo_table, geo_cache, geo_provider, veg_indices, compute_std,
                      outfile, keep_going, progress, verbose, command) where {T}
    output_options = _output_options(config.options)
    isempty(output_options.statistics) ||
        throw(ArgumentError("output.statistics requires L2 footprint gridding"))
    metadata = InputMetadata(config)
    dates = collect(time_spec.start_date:time_spec.time_step:time_spec.stop_date)
    n_times = length(dates)
    n_lon = length(grid_spec.lon)
    n_lat = length(grid_spec.lat)
    n_vars = length(config.grid_vars)

    compute_std && error("grid_center does not support compute_std yet")
    is_rolling(time_spec) && error("Rolling-mean mode is not supported for center gridding yet; use the l2 command, or drop --sample_step.")

    # Extra variables for vegetation indices
    veg_names = veg_indices ? _VEGETATION_NAMES : ()
    n_total = n_vars + length(veg_names)

    verbose && @info "Output file dimension (time/lon/lat): $n_times/$n_lon/$n_lat"

    provider = center_geolocation_provider(config;
                                           geo_table=geo_table,
                                           geo_cache=geo_cache,
                                           geo_provider=geo_provider)

    ds_out = nothing

    total_files = 0
    successful_files = Ref(0)
    failed_files = 0

    try
        open_geolocation!(provider)
        ds_out, nc_vars = create_output_dataset(outfile, grid_spec, n_times,
                                                config.grid_vars, compute_std;
                                                command=command,
                                                filters=config.filters,
                                                output_options=output_options,
                                                vegetation_names=veg_names)

        grid_data = zeros(T, n_lon, n_lat, n_total)
        grid_weights = zeros(T, n_lon, n_lat)

        p_time = _progress_bar(n_times, "Time steps: ", progress)

        # Per-tile cache: when consecutive files share a tile (file sort puts
        # same-tile granules adjacent), `lat_in`/`lon_in` come back as the
        # exact same arrays from the geolocation provider's single-slot cache.
        # We can then skip recomputing `idx`, `lat_idx`, `lon_idx`, which is
        # the same ~150 ms of bbox filter and grid-index math for every tile
        # repeat. Identity (===) check works because only sinusoidal-MODIS
        # caches; basic providers always return fresh arrays.
        cache = _center_cache(provider)
        work = (; config, grid_spec, provider, cache, veg_indices, nc_vars, metadata,
                copy_attributes=Ref(true), successful_files)

        for (ct, d) in enumerate(dates)
            ProgressMeter.next!(p_time; showvalues=[(:Time, d)])

            end_date = d + time_spec.time_step - Dates.Day(1)
            files = find_files_range(config.file_pattern, config.folder, d, end_date)
            files = sort_files_for_provider(provider, files)
            total_files += length(files)

            p_files = _progress_bar(length(files), "  Files: ", progress)

            for filepath in files
                try
                    _process_center_file!(grid_data, grid_weights, work, filepath, p_files)
                catch e
                    failed_files += 1
                    if keep_going
                        @warn "Error processing file; continuing" filepath error=sprint(showerror, e)
                    else
                        error("Error processing file $filepath: $(sprint(showerror, e))")
                    end
                end
            end

            # Write output for this time step
            min_count = T(get(config.options, "min_count", 5))
            _write_center_time_slice!(ds_out, nc_vars, config.grid_vars, veg_names, output_options,
                                      grid_data, grid_weights, ct, d, min_count)

            fill!(grid_data, zero(T))
            fill!(grid_weights, zero(T))
        end

        total_files > 0 || error("No input files matched pattern '$(config.file_pattern)' in folder '$(config.folder)'")
        successful_files[] > 0 || error("No input files were successfully processed; failed files: $failed_files")
    finally
        try
            close_geolocation!(provider)
        finally
            _close_output_dataset!(ds_out, metadata)
        end
    end

    verbose && @info "Output written to: $outfile"
    nothing
end

function _process_center_file!(grid_data, grid_weights, work, filepath, progress)
    (; config, grid_spec, provider, cache, veg_indices, nc_vars, metadata,
       copy_attributes, successful_files) = work

    # Cheap rejection of tiles that can't intersect the grid bbox. For MODIS
    # sinusoidal we know the extent analytically, so neither the per-tile cache
    # file nor the granule's HDF has to be read.
    if !provider_tile_intersects_grid(provider, filepath, grid_spec)
        ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, 0)])
        successful_files[] += 1
        return nothing
    end

    lat, lon = center_coordinates(provider, filepath, config)
    entry = _center_indices!(cache, lat, lon, grid_spec)
    if entry === nothing
        ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, 0)])
        successful_files[] += 1
        return nothing
    end

    values, valid = _center_values(filepath, config, entry.idx, veg_indices, eltype(grid_data))
    n_valid = count(valid)
    if n_valid > 0
        accumulate_center!(grid_data, grid_weights, entry.lat_idx, entry.lon_idx,
                           values, valid, length(entry.lat_idx), size(grid_data, 3))
    end
    if !_looks_like_hdf4(filepath)
        Dataset(filepath) do fin
            _capture_input_metadata!(metadata, fin)
            if copy_attributes[]
                copy_variable_attributes!(nc_vars, fin, config.grid_vars;
                    attributes=get(metadata.settings, "variable_attributes", DEFAULT_VARIABLE_ATTRIBUTES))
                copy_attributes[] = false
            end
        end
    end
    successful_files[] += 1
    ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, length(entry.idx)),
                                              (:N_valid, n_valid)])
    nothing
end

function _center_indices!(cache, lat, lon, spec)
    entry = cache[]
    if entry !== nothing && lat === entry.lat && lon === entry.lon
        # A cached empty selection still reaches input and metadata handling.
        return entry
    end
    idx = apply_center_filters(lat, lon, spec)
    lat_idx, lon_idx = _center_grid_indices(idx, lat, lon, spec)
    entry = _CenterIndexCache(lat, lon, idx, lat_idx, lon_idx)
    cache[] = entry
    isempty(idx) ? nothing : entry
end

function _center_grid_indices(idx, lat, lon, grid_spec::RectangularGridSpec)
    lat_idx = Vector{Int32}(undef, length(idx))
    lon_idx = Vector{Int32}(undef, length(idx))
    n_lat = Int32(length(grid_spec.lat))
    n_lon = Int32(length(grid_spec.lon))
    @inbounds for (i, ci) in enumerate(idx)
        ilat = floor(Int32, (lat[ci] - grid_spec.lat_min) / grid_spec.dlat) + Int32(1)
        ilon = floor(Int32, (lon[ci] - grid_spec.lon_min) / grid_spec.dlon) + Int32(1)
        lat_idx[i] = clamp(ilat, Int32(1), n_lat)
        lon_idx[i] = clamp(ilon, Int32(1), n_lon)
    end
    lat_idx, lon_idx
end

function _center_values(filepath::String, config::DataSourceConfig,
                        idx, veg_indices::Bool, ::Type{T}) where {T}
    n_pixels = length(idx)
    n_vars = length(config.grid_vars)
    n_veg = veg_indices ? length(_VEGETATION_NAMES) : 0
    values = Matrix{T}(undef, n_pixels, n_vars + n_veg)
    valid = trues(n_pixels)

    scale::T = T(get(config.options, "scale_factor", 1.0))
    offset::T = T(get(config.options, "add_offset", 0.0))
    fill_value = get(config.options, "fill_value", nothing)
    valid_min = get(config.options, "valid_min", nothing)
    valid_max = get(config.options, "valid_max", nothing)
    transpose_data::Bool = Bool(get(config.options, "transpose_data", false))

    varpaths = [varpath for (_, varpath) in config.grid_vars]
    arrays = read_arrays_from_file(filepath, vcat(varpaths, [r.variable for r in config.filters]))
    for co in 1:n_vars
        # Function-barrier on `data` so the per-pixel inner loop specializes on
        # the array element type instead of dispatching dynamically per call.
        _fill_center_column!(values, valid, idx, arrays[co], co,
                             scale, offset, fill_value, valid_min, valid_max,
                             transpose_data)
    end

    for (j, rule) in enumerate(config.filters)
        data = arrays[n_vars + j]
        for (i, ci) in enumerate(idx)
            raw = _center_data_value(data, ci, transpose_data)
            valid[i] &= _passes_quality_filter(raw, rule)
        end
    end

    nir_index = _center_var_position(config, "vegetation_nir", 2;
                                     require_default=false)
    min_nir = get(config.options, "min_nir_reflectance", nothing)
    if min_nir !== nothing && n_vars >= 2
        min_nir_t::T = T(min_nir)
        @inbounds for i in 1:n_pixels
            if values[i, nir_index] <= min_nir_t
                valid[i] = false
            end
        end
    end

    if veg_indices
        red_index = _center_var_position(config, "vegetation_red", 1)
        nir_index = _center_var_position(config, "vegetation_nir", 2)
        blue_index = _center_var_position(config, "vegetation_blue", 3)
        swir_index = _center_var_position(config, "vegetation_swir", 5)
        @inbounds for i in 1:n_pixels
            red = values[i, red_index]
            nir = values[i, nir_index]
            blue = values[i, blue_index]
            swir = values[i, swir_index]
            values[i, n_vars + 1] = compute_evi(red, nir, blue; L=T(1))
            values[i, n_vars + 2] = compute_ndvi(red, nir)
            values[i, n_vars + 3] = compute_nirv(red, nir)
            values[i, n_vars + 4] = compute_ndwi(nir, swir)
        end
        _mark_finite_center_rows!(valid, values, n_pixels, n_vars + n_veg)
    end

    values, valid
end

function _center_var_position(config::DataSourceConfig, option::String,
                              default_index::Integer; require_default::Bool=true)
    grid_keys = collect(keys(config.grid_vars))
    if haskey(config.options, option)
        var_name = string(config.options[option])
        idx = findfirst(==(var_name), grid_keys)
        idx !== nothing ||
            error("Config option '$option' references '$var_name', but that key is not in [grid]")
        return idx
    end

    if length(grid_keys) >= default_index
        return Int(default_index)
    end

    require_default &&
        error("Missing config option '$option' and [grid] has fewer than $default_index variables")
    return min(length(grid_keys), Int(default_index))
end

@inline function _center_data_value(data, ci, transpose_data::Bool)
    if transpose_data && ci isa CartesianIndex{2}
        return data[ci[2], ci[1]]
    end
    data[ci]
end

# Narrow Any-typed config bounds to the data element type where possible, while
# preserving the original value if conversion would overflow or change type.
@inline _typed_bound(::Type{D}, ::Nothing) where {D} = nothing
@inline function _typed_bound(::Type{D}, x) where {D}
    try
        return convert(D, x)
    catch
        return x
    end
end

function _fill_center_column!(values::AbstractMatrix{T}, valid::BitVector,
                              idx, data::AbstractArray, co::Integer,
                              scale::T, offset::T,
                              fill_value, valid_min, valid_max,
                              transpose_data::Bool) where {T}
    D = eltype(data)
    fv = _typed_bound(D, fill_value)
    lo = _typed_bound(D, valid_min)
    hi = _typed_bound(D, valid_max)
    @inbounds for (i, ci) in enumerate(idx)
        raw = _center_data_value(data, ci, transpose_data)
        if _invalid_center_value(raw, fv, lo, hi)
            valid[i] = false
            values[i, co] = T(NaN)
        else
            values[i, co] = T(raw) * scale + offset
        end
    end
    nothing
end

function _mark_finite_center_rows!(valid::BitVector, values::AbstractMatrix,
                                   n_pixels::Int, n_cols::Int)
    @inbounds for i in 1:n_pixels
        row_ok = true
        for co in 1:n_cols
            if !isfinite(values[i, co])
                row_ok = false
                break
            end
        end
        valid[i] &= row_ok
    end
    nothing
end

@inline function _invalid_center_value(raw, fill_value, valid_min, valid_max)
    ismissing(raw) && return true
    fill_value !== nothing && raw == fill_value && return true
    valid_min !== nothing && raw < valid_min && return true
    valid_max !== nothing && raw > valid_max && return true
    return !isfinite(raw)
end
