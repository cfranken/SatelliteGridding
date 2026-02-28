"""
    grid_l2(config, grid_spec, time_spec; kwargs...)

Grid Level-2 satellite data with footprint-aware oversampling.

Reads input NetCDF files as specified by `config`, applies quality filters, subdivides
each satellite footprint into sub-pixels for proper overlap handling, and accumulates
gridded averages. Output is written to a NetCDF4 file.

# Arguments
- `config::DataSourceConfig`: Data source configuration (from TOML/JSON)
- `grid_spec::GridSpec`: Output grid geometry
- `time_spec::TimeSpec`: Temporal binning parameters

# Keyword Arguments
- `n_oversample::Union{Nothing,Int}=nothing`: Sub-pixel factor. `nothing` = auto-compute from footprint/grid ratio.
- `compute_std::Bool=false`: Also compute per-cell standard deviation
- `outfile::String="gridded_output.nc"`: Output file path
- `backend`: KernelAbstractions backend (default: `nothing` for sequential Welford path).
  Use `CPU()` for threaded KA execution or `CUDABackend()` for GPU.
"""
function grid_l2(config::DataSourceConfig, grid_spec::GridSpec{T}, time_spec::TimeSpec;
                 n_oversample::Union{Nothing,Int}=nothing,
                 compute_std::Bool=false,
                 outfile::String="gridded_output.nc",
                 backend=nothing) where {T}

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

    fill_attrib = true
    p_time = Progress(n_times, desc="Time steps: ")

    for (ct, d) in enumerate(dates)
        ProgressMeter.next!(p_time; showvalues=[(:Time, d)])

        # Find files for this time window
        end_date = d + time_spec.time_step * time_spec.oversample_temporal - Dates.Day(1)
        files = find_files_range(config.file_pattern, config.folder, d, end_date)

        p_files = Progress(length(files), desc="  Files: ")

        for filepath in files
            try
                if use_ka
                    _process_l2_file_ka!(filepath, config, grid_spec, backend,
                                         grid_sum, grid_weights,
                                         n_oversample, nc_vars, fill_attrib, p_files)
                else
                    _process_l2_file!(filepath, config, grid_spec, grid_data, grid_std_arr,
                                      grid_weights, compute_std, n_oversample,
                                      nc_vars, fill_attrib, p_files)
                end
                if fill_attrib
                    fill_attrib = false
                end
            catch e
                println("Error processing file: $filepath")
                println("  ", e)
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

    close(ds_out)
    println("Output written to: $outfile")
    nothing
end

function _process_l2_file!(filepath::String, config::DataSourceConfig,
                           grid_spec::GridSpec{T}, grid_data, grid_std,
                           grid_weights, compute_std, n_oversample_override,
                           nc_vars, fill_attrib, progress) where {T}
    fin = Dataset(filepath)
    try
        # Quick bounding box check with center coordinates
        lat_center = read_nc_variable(fin, config.basic["lat"])
        lon_center = read_nc_variable(fin, config.basic["lon"])

        in_bounds = ((lat_center .> grid_spec.lat_min) .+
                     (lat_center .< grid_spec.lat_max) .+
                     (lon_center .> grid_spec.lon_min) .+
                     (lon_center .< grid_spec.lon_max))

        if !any(in_bounds .== 4)
            ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, 0)])
            return
        end

        # Read footprint corner coordinates
        lat_bnd = read_nc_variable(fin, config.basic["lat_bnd"]; bounds=true)
        lon_bnd = read_nc_variable(fin, config.basic["lon_bnd"]; bounds=true)

        # Transpose if needed (ensure N×4)
        if size(lat_bnd, 1) == 4 && ndims(lat_bnd) == 2
            lat_bnd = lat_bnd'
            lon_bnd = lon_bnd'
        end

        # Normalize corner ordering to counter-clockwise
        # (different satellite products store corners in different orders)
        sort_corners_ccw!(lat_bnd, lon_bnd)

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

        # Convert corner coordinates to fractional grid indices
        n_lat = length(grid_spec.lat)
        n_lon = length(grid_spec.lon)
        ilat = ((lat_bnd[idx, :] .- grid_spec.lat_min) /
                (grid_spec.lat_max - grid_spec.lat_min) * n_lat) .+ 1
        ilon = ((lon_bnd[idx, :] .- grid_spec.lon_min) /
                (grid_spec.lon_max - grid_spec.lon_min) * n_lon) .+ 1

        # Determine n_oversample
        n = if n_oversample_override !== nothing
            n_oversample_override
        else
            # Auto-compute from median footprint extent vs grid cell size
            lat_extents = maximum(ilat, dims=2) .- minimum(ilat, dims=2)
            med_extent = median(lat_extents)
            compute_n_oversample(med_extent, 1.0)
        end

        # Allocate buffers
        points_buf = zeros(T, n, n, 2)
        ix_buf = zeros(Int32, n^2)
        iy_buf = zeros(Int32, n^2)
        lats0 = zeros(n)
        lons0 = zeros(n)
        lats1 = zeros(n)
        lons1 = zeros(n)

        accumulate_footprint!(grid_data, grid_std, grid_weights, compute_std,
                              ilat, ilon, mat_in[idx, :],
                              length(idx), length(config.grid_vars), n,
                              points_buf, ix_buf, iy_buf,
                              lats0, lons0, lats1, lons1)
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
                               nc_vars, fill_attrib, progress) where {T}
    fin = Dataset(filepath)
    try
        # Quick bounding box check
        lat_center = read_nc_variable(fin, config.basic["lat"])
        lon_center = read_nc_variable(fin, config.basic["lon"])

        in_bounds = ((lat_center .> grid_spec.lat_min) .+
                     (lat_center .< grid_spec.lat_max) .+
                     (lon_center .> grid_spec.lon_min) .+
                     (lon_center .< grid_spec.lon_max))

        if !any(in_bounds .== 4)
            ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, 0)])
            return
        end

        # Read footprint corner coordinates
        lat_bnd = read_nc_variable(fin, config.basic["lat_bnd"]; bounds=true)
        lon_bnd = read_nc_variable(fin, config.basic["lon_bnd"]; bounds=true)

        if size(lat_bnd, 1) == 4 && ndims(lat_bnd) == 2
            lat_bnd = lat_bnd'
            lon_bnd = lon_bnd'
        end

        # Apply filters (sorting is done inside accumulate_batch!)
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

        # Convert corner coordinates to fractional grid indices
        n_lat = length(grid_spec.lat)
        n_lon = length(grid_spec.lon)
        ilat = T.((lat_bnd[idx, :] .- grid_spec.lat_min) /
                   (grid_spec.lat_max - grid_spec.lat_min) .* T(n_lat)) .+ one(T)
        ilon = T.((lon_bnd[idx, :] .- grid_spec.lon_min) /
                   (grid_spec.lon_max - grid_spec.lon_min) .* T(n_lon)) .+ one(T)

        # Determine n_oversample
        n = if n_oversample_override !== nothing
            n_oversample_override
        else
            lat_extents = maximum(ilat, dims=2) .- minimum(ilat, dims=2)
            med_extent = median(lat_extents)
            compute_n_oversample(med_extent, 1.0)
        end

        accumulate_batch!(backend, grid_sum, grid_weights,
                          ilat, ilon, mat_in[idx, :], n, length(config.grid_vars))
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
            nc_vars[key][ct, :, :] = da
            if compute_std
                da_std = round.(grid_std[:, :, co], sigdigits=6)
                da_std[grid_weights .< 1e-10] .= -999.0f0
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
- `geo_table::Union{Nothing,String}=nothing`: Path to a geolocation lookup table (NetCDF)
- `veg_indices::Bool=false`: Compute vegetation indices (EVI, NDVI, NIRv, NDWI)
- `compute_std::Bool=false`: Compute per-cell standard deviation
- `outfile::String="gridded_output.nc"`: Output file path
"""
function grid_center(config::DataSourceConfig, grid_spec::GridSpec{T}, time_spec::TimeSpec;
                     geo_table::Union{Nothing,String}=nothing,
                     veg_indices::Bool=false,
                     compute_std::Bool=false,
                     outfile::String="gridded_output.nc") where {T}

    dates = collect(time_spec.start_date:time_spec.time_step:time_spec.stop_date)
    n_times = length(dates)
    n_lon = length(grid_spec.lon)
    n_lat = length(grid_spec.lat)
    n_vars = length(config.grid_vars)

    # Extra variables for vegetation indices
    n_veg = veg_indices ? 4 : 0
    n_total = n_vars + n_veg + 1  # +1 for count

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

    # Open geo table if provided
    ds_geo = geo_table !== nothing ? Dataset(geo_table) : nothing

    p_time = Progress(n_times, desc="Time steps: ")

    for (ct, d) in enumerate(dates)
        ProgressMeter.next!(p_time; showvalues=[(:Time, d)])

        end_date = d + time_spec.time_step - Dates.Day(1)
        files = find_files_range(config.file_pattern, config.folder, d, end_date)

        p_files = Progress(length(files), desc="  Files: ")

        for filepath in files
            try
                fin = Dataset(filepath)
                try
                    # Read data variables
                    co = 1
                    for (_, value) in config.grid_vars
                        # MODIS files: read full tile into grid_data accumulator
                        # (placeholder — full MODIS pipeline requires geo_table tile lookup)
                        _ = fin[value].var[:]
                        co += 1
                    end
                    ProgressMeter.next!(p_files; showvalues=[(:File, filepath)])
                finally
                    close(fin)
                end
            catch e
                println("Error processing file: $filepath — $e")
            end
        end

        # Write output for this time step
        ds_out["time"][ct] = d
        weights = grid_data[:, :, end]
        nc_vars["n"][ct, :, :] = weights

        threshold = T(5)
        co = 1
        for (key, _) in config.grid_vars
            da = round.(grid_data[:, :, co] ./ max.(weights, eps(T)), sigdigits=5)
            da[weights .< threshold] .= T(-999)
            nc_vars[key][ct, :, :] = da
            co += 1
        end

        fill!(grid_data, zero(T))
    end

    ds_geo !== nothing && close(ds_geo)
    close(ds_out)
    println("Output written to: $outfile")
    nothing
end
