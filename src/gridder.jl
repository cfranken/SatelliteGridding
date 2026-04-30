"""
    grid_l2(config, grid_spec, time_spec; kwargs...)

Grid Level-2 satellite data with footprint-aware oversampling.

Reads input NetCDF files as specified by `config`, applies quality filters, subdivides
each satellite footprint into sub-pixels for proper overlap handling, and accumulates
gridded averages. Output is written to a NetCDF4 file. Quadrilateral footprints
are used by default; pass `footprint_method=CircularFootprintGridding(...)` for
circular products such as GOSAT. Circular footprints can be read from
`lat_bnd`/`lon_bnd` or from center `lat`/`lon` plus `radius`.

# Arguments
- `config::DataSourceConfig`: Data source configuration (from TOML/JSON)
- `grid_spec::GridSpec`: Output grid geometry
- `time_spec::TimeSpec`: Temporal binning parameters

# Keyword Arguments
- `n_oversample::Union{Nothing,Int}=nothing`: Sub-pixel factor. `nothing` = auto-compute from footprint/grid ratio.
- `footprint_method::AbstractGriddingMethod=SubpixelGridding(n_oversample)`: Footprint geometry method.
- `compute_std::Bool=false`: Also compute per-cell standard deviation
- `outfile::String="gridded_output.nc"`: Output file path
- `backend`: KernelAbstractions backend (default: `nothing` for sequential Welford path).
  Use `CPU()`, `CUDABackend()`, `MetalBackend()`, or another compatible KA backend.
- `keep_going::Bool=false`: Continue after per-file errors. By default, the
  first file error stops the run.
"""
function grid_l2(config::DataSourceConfig, grid_spec::GridSpec{T}, time_spec::TimeSpec;
                 n_oversample::Union{Nothing,Int}=nothing,
                 footprint_method::AbstractGriddingMethod=SubpixelGridding(n_oversample),
                 compute_std::Bool=false,
                 outfile::String="gridded_output.nc",
                 backend=nothing,
                 keep_going::Bool=false) where {T}

    footprint_method isa SubpixelGridding ||
        footprint_method isa CircularFootprintGridding ||
        error("grid_l2 supports SubpixelGridding and CircularFootprintGridding; use grid_center for center-coordinate gridding")
    _validate_l2_geometry_config(config, footprint_method)

    effective_n_oversample = if footprint_method isa Union{SubpixelGridding,CircularFootprintGridding} &&
                                footprint_method.n_oversample !== nothing
        footprint_method.n_oversample
    else
        n_oversample
    end

    dates = collect(time_spec.start_date:time_spec.time_step:time_spec.stop_date)
    n_times = length(dates)
    n_lon = length(grid_spec.lon)
    n_lat = length(grid_spec.lat)
    n_vars = length(config.grid_vars)

    use_ka = backend !== nothing

    println("Output file dimension (time/lon/lat): $n_times/$n_lon/$n_lat")
    use_ka && println("Using KA backend: $backend")

    # Create output file
    ds_out, nc_vars = create_output_dataset(outfile, grid_spec, n_times,
                                            config.grid_vars, compute_std)

    # Allocate work arrays (on backend for GPU, regular Array for CPU)
    if use_ka
        grid_sum = KernelAbstractions.zeros(backend, T, n_lon, n_lat, n_vars)
        grid_weights = KernelAbstractions.zeros(backend, T, n_lon, n_lat)
    else
        # Welford running mean for sequential path
        grid_data = zeros(T, n_lon, n_lat, n_vars)
        grid_std_arr = zeros(T, n_lon, n_lat, n_vars)
        grid_weights = zeros(T, n_lon, n_lat)
    end

    total_files = 0
    successful_files = 0
    failed_files = 0

    try
        fill_attrib = true
        p_time = Progress(n_times, desc="Time steps: ")

        for (ct, d) in enumerate(dates)
            ProgressMeter.next!(p_time; showvalues=[(:Time, d)])

            # Find files for this time window
            end_date = d + time_spec.time_step * time_spec.oversample_temporal - Dates.Day(1)
            files = find_files_range(config.file_pattern, config.folder, d, end_date)
            total_files += length(files)

            p_files = Progress(length(files), desc="  Files: ")

            for filepath in files
                try
                    if use_ka
                        _process_l2_file_ka!(filepath, config, grid_spec, backend,
                                             grid_sum, grid_weights,
                                             effective_n_oversample,
                                             footprint_method,
                                             nc_vars, fill_attrib, p_files)
                    else
                        _process_l2_file!(filepath, config, grid_spec, grid_data, grid_std_arr,
                                          grid_weights, compute_std,
                                          effective_n_oversample,
                                          footprint_method,
                                          nc_vars, fill_attrib, p_files)
                    end
                    successful_files += 1
                    if fill_attrib
                        fill_attrib = false
                    end
                catch e
                    failed_files += 1
                    if keep_going
                        @warn "Error processing file; continuing" filepath error=sprint(showerror, e)
                    else
                        error("Error processing file $filepath: $(sprint(showerror, e))")
                    end
                end
            end

            # Write time slice to output
            if use_ka
                # Copy to CPU for finalization and NetCDF write (no-op for CPU backend)
                grid_data_final = Array(grid_sum)
                grid_w_cpu = Array(grid_weights)
                finalize_mean!(grid_data_final, grid_w_cpu)
                _write_time_slice!(ds_out, nc_vars, config, grid_data_final,
                                   zeros(T, n_lon, n_lat, n_vars),
                                   grid_w_cpu, false, ct, d)
                fill!(grid_sum, zero(T))
            else
                if compute_std
                    finalize_std!(grid_std_arr, grid_weights)
                end
                _write_time_slice!(ds_out, nc_vars, config, grid_data, grid_std_arr,
                                   grid_weights, compute_std, ct, d)
                fill!(grid_data, zero(T))
                fill!(grid_std_arr, zero(T))
            end

            # Reset accumulators
            fill!(grid_weights, zero(T))
        end

        total_files > 0 || error("No input files matched pattern '$(config.file_pattern)' in folder '$(config.folder)'")
        successful_files > 0 || error("No input files were successfully processed; failed files: $failed_files")
    finally
        close(ds_out)
    end

    println("Output written to: $outfile")
    nothing
end

_has_l2_bounds(config::DataSourceConfig) =
    haskey(config.basic, "lat_bnd") && haskey(config.basic, "lon_bnd")

_has_l2_centers(config::DataSourceConfig) =
    haskey(config.basic, "lat") && haskey(config.basic, "lon")

_has_l2_radius(config::DataSourceConfig) =
    haskey(config.basic, "radius") || haskey(config.options, "radius")

function _use_radius_l2_geometry(config::DataSourceConfig,
                                 footprint_method::AbstractGriddingMethod)
    footprint_method isa CircularFootprintGridding &&
        _has_l2_radius(config) && !_has_l2_bounds(config)
end

function _validate_l2_geometry_config(config::DataSourceConfig,
                                      footprint_method::AbstractGriddingMethod)
    has_bounds = _has_l2_bounds(config)
    has_centers = _has_l2_centers(config)
    has_radius = _has_l2_radius(config)

    has_lat = haskey(config.basic, "lat")
    has_lon = haskey(config.basic, "lon")
    has_lat == has_lon ||
        error("grid_l2 requires both basic.lat and basic.lon, or neither so centers can be derived from bounds")

    if footprint_method isa CircularFootprintGridding
        has_bounds || has_radius ||
            error("CircularFootprintGridding requires either basic.lat_bnd/basic.lon_bnd or center basic.lat/basic.lon plus basic.radius or [circle] radius")
        if !has_bounds
            has_centers ||
                error("Radius-based CircularFootprintGridding requires basic.lat and basic.lon center variables")
        elseif !has_centers
            @warn "basic.lat/basic.lon not configured; deriving footprint centers from corner coordinates"
        end
    else
        has_bounds ||
            error("grid_l2 requires basic.lat_bnd and basic.lon_bnd for quadrilateral footprint gridding")
        if !has_centers
            @warn "basic.lat/basic.lon not configured; deriving footprint centers from corner coordinates"
        end
    end
    nothing
end

function _read_or_derive_l2_centers(fin, config::DataSourceConfig,
                                    lat_bnd::AbstractMatrix,
                                    lon_bnd::AbstractMatrix)
    if haskey(config.basic, "lat") && haskey(config.basic, "lon")
        return vec(read_nc_variable(fin, config.basic["lat"])),
               vec(read_nc_variable(fin, config.basic["lon"]))
    end

    vec(mean(lat_bnd, dims=2)), vec(mean(lon_bnd, dims=2))
end

function _read_l2_centers(fin, config::DataSourceConfig)
    _has_l2_centers(config) ||
        error("Radius-based circular footprint gridding requires basic.lat and basic.lon")
    vec(read_nc_variable(fin, config.basic["lat"])),
    vec(read_nc_variable(fin, config.basic["lon"]))
end

function _read_l2_radius(fin, config::DataSourceConfig, n::Integer)
    if haskey(config.basic, "radius")
        radius = vec(read_nc_variable(fin, config.basic["radius"]))
        length(radius) == n ||
            error("Radius variable '$(config.basic["radius"])' has $(length(radius)) values, expected $n")
        return radius
    end

    haskey(config.options, "radius") ||
        error("Radius-based circular footprint gridding requires basic.radius or [circle] radius")
    fill(Float64(config.options["radius"]), n)
end

function _radius_to_degrees(radius::Real, lat::Real, unit::AbstractString)
    unit_lc = lowercase(strip(unit))
    if unit_lc in ("degree", "degrees", "deg")
        r = Float64(radius)
    elseif unit_lc in ("kilometer", "kilometers", "km")
        r = Float64(radius) / 111.32
    elseif unit_lc in ("meter", "meters", "m")
        r = Float64(radius) / 111_320.0
    else
        error("Unknown [circle] radius_unit '$unit'. Use degrees, km, or m.")
    end

    lon_scale = max(abs(cosd(Float64(lat))), 1e-6)
    r, r / lon_scale
end

function _circle_bounds_from_radius(lat_center, lon_center, radius,
                                    config::DataSourceConfig)
    unit = string(get(config.options, "radius_unit", "degrees"))
    n = length(lat_center)
    lat_bnd = zeros(Float64, n, 4)
    lon_bnd = zeros(Float64, n, 4)

    @inbounds for i in 1:n
        if !isfinite(radius[i]) || radius[i] < 0
            lat_bnd[i, :] .= NaN
            lon_bnd[i, :] .= NaN
            continue
        end
        rlat, rlon = _radius_to_degrees(radius[i], lat_center[i], unit)
        lat_bnd[i, 1] = lat_center[i] - rlat
        lat_bnd[i, 2] = lat_center[i] - rlat
        lat_bnd[i, 3] = lat_center[i] + rlat
        lat_bnd[i, 4] = lat_center[i] + rlat
        lon_bnd[i, 1] = lon_center[i] - rlon
        lon_bnd[i, 2] = lon_center[i] + rlon
        lon_bnd[i, 3] = lon_center[i] + rlon
        lon_bnd[i, 4] = lon_center[i] - rlon
    end

    lat_bnd, lon_bnd
end

function _read_l2_geometry(fin, config::DataSourceConfig,
                           footprint_method::AbstractGriddingMethod)
    if _use_radius_l2_geometry(config, footprint_method)
        lat_center, lon_center = _read_l2_centers(fin, config)
        radius = _read_l2_radius(fin, config, length(lat_center))
        lat_bnd, lon_bnd = _circle_bounds_from_radius(lat_center, lon_center,
                                                      radius, config)
        return lat_center, lon_center, lat_bnd, lon_bnd
    end

    lat_bnd = read_nc_variable(fin, config.basic["lat_bnd"]; bounds=true)
    lon_bnd = read_nc_variable(fin, config.basic["lon_bnd"]; bounds=true)

    if size(lat_bnd, 1) == 4 && ndims(lat_bnd) == 2
        lat_bnd = lat_bnd'
        lon_bnd = lon_bnd'
    end

    sort_corners_ccw!(lat_bnd, lon_bnd)
    lat_center, lon_center = _read_or_derive_l2_centers(fin, config,
                                                        lat_bnd, lon_bnd)
    lat_center, lon_center, lat_bnd, lon_bnd
end

function _process_l2_file!(filepath::String, config::DataSourceConfig,
                           grid_spec::GridSpec{T}, grid_data, grid_std,
                           grid_weights, compute_std, n_oversample_override,
                           footprint_method::AbstractGriddingMethod,
                           nc_vars, fill_attrib, progress) where {T}
    fin = Dataset(filepath)
    try
        lat_center, lon_center, lat_bnd, lon_bnd =
            _read_l2_geometry(fin, config, footprint_method)

        # Quick bounding box check with center coordinates.
        in_bounds = ((lat_center .> grid_spec.lat_min) .+
                     (lat_center .< grid_spec.lat_max) .+
                     (lon_center .> grid_spec.lon_min) .+
                     (lon_center .< grid_spec.lon_max))

        if !any(in_bounds .== 4)
            ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, 0)])
            return
        end

        # Apply filters
        idx = apply_filters(fin, config, lat_bnd, lon_bnd, grid_spec)
        ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, length(idx))])

        if isempty(idx)
            return
        end

        # Copy attributes from first file
        if fill_attrib
            copy_variable_attributes!(nc_vars, fin, config.grid_vars)
        end

        # Read all grid variables
        n_soundings = size(lat_bnd, 1)
        mat_in = zeros(T, n_soundings, length(config.grid_vars))
        co = 1
        for (_, value) in config.grid_vars
            mat_in[:, co] = read_nc_variable(fin, value)
            co += 1
        end

        # Extract values and corners for filtered soundings
        vals = mat_in[idx, :]
        lat_c = lat_bnd[idx, :]
        lon_c = lon_bnd[idx, :]
        lat_center_c = lat_center[idx]
        lon_center_c = lon_center[idx]

        # Filter out soundings with NaN/Inf in data or corners
        finite_mask = vec(all(isfinite, vals, dims=2)) .&
                      vec(all(isfinite, lat_c, dims=2)) .&
                      vec(all(isfinite, lon_c, dims=2)) .&
                      isfinite.(lat_center_c) .&
                      isfinite.(lon_center_c)
        if !all(finite_mask)
            vals = vals[finite_mask, :]
            lat_c = lat_c[finite_mask, :]
            lon_c = lon_c[finite_mask, :]
            lat_center_c = lat_center_c[finite_mask]
            lon_center_c = lon_center_c[finite_mask]
        end

        if isempty(vals)
            return
        end

        # Convert corner coordinates to fractional grid indices
        n_lat = length(grid_spec.lat)
        n_lon = length(grid_spec.lon)
        ilat = ((lat_c .- grid_spec.lat_min) /
                (grid_spec.lat_max - grid_spec.lat_min) * n_lat) .+ 1
        ilon = ((lon_c .- grid_spec.lon_min) /
                (grid_spec.lon_max - grid_spec.lon_min) * n_lon) .+ 1
        ilat_center = ((lat_center_c .- grid_spec.lat_min) /
                       (grid_spec.lat_max - grid_spec.lat_min) * n_lat) .+ 1
        ilon_center = ((lon_center_c .- grid_spec.lon_min) /
                       (grid_spec.lon_max - grid_spec.lon_min) * n_lon) .+ 1

        # Determine n_oversample
        n = if n_oversample_override !== nothing
            n_oversample_override
        else
            # Auto-compute from median footprint extent vs grid cell size
            lat_extents = maximum(ilat, dims=2) .- minimum(ilat, dims=2)
            lon_extents = maximum(ilon, dims=2) .- minimum(ilon, dims=2)
            med_extent = median(max.(lat_extents, lon_extents))
            compute_n_oversample(med_extent, 1.0)
        end

        if footprint_method isa CircularFootprintGridding
            accumulate_circular_footprint!(grid_data, grid_std, grid_weights,
                                           compute_std, ilat_center, ilon_center,
                                           ilat, ilon, vals,
                                           size(vals, 1), length(config.grid_vars), n)
        else
            # Allocate buffers
            points_buf = zeros(T, n, n, 2)
            ix_buf = zeros(Int32, n^2)
            iy_buf = zeros(Int32, n^2)
            lats0 = zeros(n)
            lons0 = zeros(n)
            lats1 = zeros(n)
            lons1 = zeros(n)

            accumulate_footprint!(grid_data, grid_std, grid_weights, compute_std,
                                  ilat, ilon, vals,
                                  size(vals, 1), length(config.grid_vars), n,
                                  points_buf, ix_buf, iy_buf,
                                  lats0, lons0, lats1, lons1)
        end
    finally
        close(fin)
    end
    nothing
end

"""
Process a single L2 file using the KA kernel pipeline (sum-based accumulation).
"""
function _process_l2_file_ka!(filepath::String, config::DataSourceConfig,
                               grid_spec::GridSpec{T}, backend,
                               grid_sum, grid_weights,
                               n_oversample_override,
                               footprint_method::AbstractGriddingMethod,
                               nc_vars, fill_attrib, progress) where {T}
    fin = Dataset(filepath)
    try
        lat_center, lon_center, lat_bnd, lon_bnd =
            _read_l2_geometry(fin, config, footprint_method)

        # Quick bounding box check with center coordinates.
        in_bounds = ((lat_center .> grid_spec.lat_min) .+
                     (lat_center .< grid_spec.lat_max) .+
                     (lon_center .> grid_spec.lon_min) .+
                     (lon_center .< grid_spec.lon_max))

        if !any(in_bounds .== 4)
            ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, 0)])
            return
        end

        # Apply filters
        idx = apply_filters(fin, config, lat_bnd, lon_bnd, grid_spec)
        ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, length(idx))])

        if isempty(idx)
            return
        end

        if fill_attrib
            copy_variable_attributes!(nc_vars, fin, config.grid_vars)
        end

        # Read all grid variables
        n_soundings = size(lat_bnd, 1)
        mat_in = zeros(T, n_soundings, length(config.grid_vars))
        co = 1
        for (_, value) in config.grid_vars
            mat_in[:, co] = read_nc_variable(fin, value)
            co += 1
        end

        # Extract values and corners for filtered soundings
        vals = mat_in[idx, :]
        lat_c = lat_bnd[idx, :]
        lon_c = lon_bnd[idx, :]
        lat_center_c = lat_center[idx]
        lon_center_c = lon_center[idx]

        # Filter out soundings with NaN/Inf in data or corners
        finite_mask = vec(all(isfinite, vals, dims=2)) .&
                      vec(all(isfinite, lat_c, dims=2)) .&
                      vec(all(isfinite, lon_c, dims=2)) .&
                      isfinite.(lat_center_c) .&
                      isfinite.(lon_center_c)
        if !all(finite_mask)
            vals = vals[finite_mask, :]
            lat_c = lat_c[finite_mask, :]
            lon_c = lon_c[finite_mask, :]
            lat_center_c = lat_center_c[finite_mask]
            lon_center_c = lon_center_c[finite_mask]
        end

        if isempty(vals)
            return
        end

        # Convert corner coordinates to fractional grid indices
        n_lat = length(grid_spec.lat)
        n_lon = length(grid_spec.lon)
        ilat = T.((lat_c .- grid_spec.lat_min) /
                   (grid_spec.lat_max - grid_spec.lat_min) .* T(n_lat)) .+ one(T)
        ilon = T.((lon_c .- grid_spec.lon_min) /
                   (grid_spec.lon_max - grid_spec.lon_min) .* T(n_lon)) .+ one(T)
        ilat_center = T.((lat_center_c .- grid_spec.lat_min) /
                          (grid_spec.lat_max - grid_spec.lat_min) .* T(n_lat)) .+ one(T)
        ilon_center = T.((lon_center_c .- grid_spec.lon_min) /
                          (grid_spec.lon_max - grid_spec.lon_min) .* T(n_lon)) .+ one(T)

        # Determine n_oversample
        n = if n_oversample_override !== nothing
            n_oversample_override
        else
            lat_extents = maximum(ilat, dims=2) .- minimum(ilat, dims=2)
            lon_extents = maximum(ilon, dims=2) .- minimum(ilon, dims=2)
            med_extent = median(max.(lat_extents, lon_extents))
            compute_n_oversample(med_extent, 1.0)
        end

        if footprint_method isa CircularFootprintGridding
            accumulate_circular_batch!(backend, grid_sum, grid_weights,
                                       ilat_center, ilon_center, ilat, ilon,
                                       vals, n, length(config.grid_vars))
        else
            accumulate_batch!(backend, grid_sum, grid_weights,
                              ilat, ilon, vals, n, length(config.grid_vars))
        end
    finally
        close(fin)
    end
    nothing
end

function _write_time_slice!(ds_out, nc_vars, config, grid_data, grid_std,
                            grid_weights, compute_std, ct, d)
    nc_vars["n"][ct, :, :] = grid_weights
    ds_out["time"][ct] = d

    if maximum(grid_weights) > 0
        co = 1
        for (key, _) in config.grid_vars
            da = round.(grid_data[:, :, co], sigdigits=8)
            da[grid_weights .< 1e-10] .= -999.0f0
            da[.!isfinite.(da)] .= -999.0f0
            nc_vars[key][ct, :, :] = da
            if compute_std
                da_std = round.(grid_std[:, :, co], sigdigits=6)
                da_std[grid_weights .< 1e-10] .= -999.0f0
                da_std[.!isfinite.(da_std)] .= -999.0f0
                nc_vars[key * "_std"][ct, :, :] = da_std
            end
            co += 1
        end
    else
        nc_vars["n"][ct, :, :] .= 0
    end
    nothing
end

"""
    grid_center(config, grid_spec, time_spec; kwargs...)

Grid data using center coordinates only (no footprint bounds), suitable for
MODIS-style data where each pixel maps to exactly one grid cell.

# Keyword Arguments
- `geo_table::Union{Nothing,String}=nothing`: Path to a legacy monolithic geolocation lookup table (NetCDF)
- `geo_cache::Union{Nothing,String}=nothing`: Directory for generated per-tile MODIS geolocation cache
- `geo_provider::Union{Symbol,String}=:auto`: `:auto`, `:variables`, `:lut`, or `:modis`
- `veg_indices::Bool=false`: Compute vegetation indices (EVI, NDVI, NIRv, NDWI)
- `compute_std::Bool=false`: Compute per-cell standard deviation
- `outfile::String="gridded_output.nc"`: Output file path
- `keep_going::Bool=false`: Continue after per-file errors. By default, the
  first file error stops the run.
"""
function grid_center(config::DataSourceConfig, grid_spec::GridSpec{T}, time_spec::TimeSpec;
                     geo_table::Union{Nothing,String}=nothing,
                     geo_cache::Union{Nothing,String}=nothing,
                     geo_provider::Union{Symbol,String}=:auto,
                     veg_indices::Bool=false,
                     compute_std::Bool=false,
                     outfile::String="gridded_output.nc",
                     keep_going::Bool=false) where {T}

    dates = collect(time_spec.start_date:time_spec.time_step:time_spec.stop_date)
    n_times = length(dates)
    n_lon = length(grid_spec.lon)
    n_lat = length(grid_spec.lat)
    n_vars = length(config.grid_vars)

    compute_std && error("grid_center does not support compute_std yet")

    # Extra variables for vegetation indices
    n_veg = veg_indices ? 4 : 0
    n_total = n_vars + n_veg

    println("Output file dimension (time/lon/lat): $n_times/$n_lon/$n_lat")

    ds_out, nc_vars = create_output_dataset(outfile, grid_spec, n_times,
                                            config.grid_vars, compute_std)

    # Add vegetation index variables if requested
    if veg_indices
        nc_vars["EVI"] = defVar(ds_out, "EVI", Float32, ("time", "lon", "lat"),
                                deflatelevel=4, fillvalue=-999.0f0)
        nc_vars["NDVI"] = defVar(ds_out, "NDVI", Float32, ("time", "lon", "lat"),
                                 deflatelevel=4, fillvalue=-999.0f0)
        nc_vars["NIRv"] = defVar(ds_out, "NIRv", Float32, ("time", "lon", "lat"),
                                 deflatelevel=4, fillvalue=-999.0f0)
        nc_vars["NDWI"] = defVar(ds_out, "NDWI", Float32, ("time", "lon", "lat"),
                                 deflatelevel=4, fillvalue=-999.0f0)
    end

    grid_data = zeros(T, n_lon, n_lat, n_total)
    grid_weights = zeros(T, n_lon, n_lat)

    provider = center_geolocation_provider(config;
                                           geo_table=geo_table,
                                           geo_cache=geo_cache,
                                           geo_provider=geo_provider)

    total_files = 0
    successful_files = 0
    failed_files = 0

    try
        open_geolocation!(provider)
        p_time = Progress(n_times, desc="Time steps: ")

        for (ct, d) in enumerate(dates)
            ProgressMeter.next!(p_time; showvalues=[(:Time, d)])

            end_date = d + time_spec.time_step - Dates.Day(1)
            files = find_files_range(config.file_pattern, config.folder, d, end_date)
            total_files += length(files)

            p_files = Progress(length(files), desc="  Files: ")

            for filepath in files
                try
                    lat_in, lon_in = center_coordinates(provider, filepath, config)
                    idx = apply_center_filters(lat_in, lon_in, grid_spec)
                    if isempty(idx)
                        ProgressMeter.next!(p_files; showvalues=[(:File, filepath), (:N_pixels, 0)])
                        successful_files += 1
                        continue
                    end

                    lat_idx, lon_idx = _center_grid_indices(idx, lat_in, lon_in, grid_spec)
                    values, valid = _center_values(filepath, config, idx, veg_indices, T)

                    if !all(valid)
                        lat_idx = lat_idx[valid]
                        lon_idx = lon_idx[valid]
                        values = values[valid, :]
                    end

                    n_valid = length(lat_idx)
                    if n_valid > 0
                        accumulate_center!(grid_data, grid_weights, lat_idx, lon_idx,
                                           values, n_valid, n_total)
                    end
                    successful_files += 1
                    ProgressMeter.next!(p_files; showvalues=[(:File, filepath),
                                                              (:N_pixels, length(idx)),
                                                              (:N_valid, n_valid)])
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
            _write_center_time_slice!(ds_out, nc_vars, config, grid_data,
                                      grid_weights, veg_indices, ct, d, min_count)

            fill!(grid_data, zero(T))
            fill!(grid_weights, zero(T))
        end

        total_files > 0 || error("No input files matched pattern '$(config.file_pattern)' in folder '$(config.folder)'")
        successful_files > 0 || error("No input files were successfully processed; failed files: $failed_files")
    finally
        close_geolocation!(provider)
        close(ds_out)
    end

    println("Output written to: $outfile")
    nothing
end

function _center_grid_indices(idx, lat, lon, grid_spec::GridSpec)
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
    n_veg = veg_indices ? 4 : 0
    values = zeros(T, n_pixels, n_vars + n_veg)
    valid = trues(n_pixels)

    scale = T(get(config.options, "scale_factor", 1.0))
    offset = T(get(config.options, "add_offset", 0.0))
    fill_value = get(config.options, "fill_value", nothing)
    valid_min = get(config.options, "valid_min", nothing)
    valid_max = get(config.options, "valid_max", nothing)
    transpose_data = Bool(get(config.options, "transpose_data", false))

    for (co, (_, varpath)) in enumerate(config.grid_vars)
        data = read_array_from_file(filepath, varpath)
        @inbounds for (i, ci) in enumerate(idx)
            raw = _center_data_value(data, ci, transpose_data)
            if _invalid_center_value(raw, fill_value, valid_min, valid_max)
                valid[i] = false
                values[i, co] = T(NaN)
            else
                values[i, co] = T(raw) * scale + offset
            end
        end
    end

    nir_index = _center_var_position(config, "vegetation_nir", 2;
                                     require_default=false)
    min_nir = get(config.options, "min_nir_reflectance", nothing)
    if min_nir !== nothing && n_vars >= 2
        min_nir_t = T(min_nir)
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
        valid .&= vec(all(isfinite, values, dims=2))
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

@inline function _invalid_center_value(raw, fill_value, valid_min, valid_max)
    ismissing(raw) && return true
    fill_value !== nothing && raw == fill_value && return true
    valid_min !== nothing && raw < valid_min && return true
    valid_max !== nothing && raw > valid_max && return true
    return !isfinite(raw)
end

function _write_center_time_slice!(ds_out, nc_vars, config, grid_data,
                                   grid_weights, veg_indices::Bool,
                                   ct::Int, d::DateTime, min_count)
    ds_out["time"][ct] = d
    nc_vars["n"][ct, :, :] = grid_weights

    co = 1
    for (key, _) in config.grid_vars
        _write_center_var!(nc_vars[key], grid_data[:, :, co],
                           grid_weights, ct, min_count)
        co += 1
    end

    if veg_indices
        for key in ("EVI", "NDVI", "NIRv", "NDWI")
            _write_center_var!(nc_vars[key], grid_data[:, :, co],
                               grid_weights, ct, min_count)
            co += 1
        end
    end
    nothing
end

function _write_center_var!(nc_var, data, weights, ct::Int, min_count)
    da = round.(data ./ max.(weights, eps(eltype(data))), sigdigits=5)
    da[weights .< min_count] .= eltype(data)(-999)
    da[.!isfinite.(da)] .= eltype(data)(-999)
    nc_var[ct, :, :] = da
    nothing
end
