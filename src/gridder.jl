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
function grid_l2(config::DataSourceConfig, grid_spec::AbstractGridSpec{T}, time_spec::TimeSpec;
                 n_oversample::Union{Nothing,Int}=nothing,
                 footprint_method::AbstractGriddingMethod=SubpixelGridding(n_oversample),
                 compute_std::Bool=false,
                 outfile::String="gridded_output.nc",
                 zarr_out::Union{Nothing,String}=nothing,
                 zarr_meta::Union{Nothing,String}=nothing,
                 zarr_epoch::Date=ZARR_EPOCH,
                 backend=nothing,
                 keep_going::Bool=false,
                 command::Union{Nothing,String}=nothing) where {T}

    footprint_method isa SubpixelGridding ||
        footprint_method isa CircularFootprintGridding ||
        error("grid_l2 supports SubpixelGridding and CircularFootprintGridding; use grid_center for center-coordinate gridding")
    _validate_l2_geometry_config(config, footprint_method)

    if backend !== nothing && grid_spec isa CubedSphereGridSpec
        error("Cubed-sphere gridding currently supports only the sequential CPU path; omit `backend` or use --backend sequential.")
    end
    if grid_spec isa CubedSphereGridSpec && footprint_method isa CircularFootprintGridding
        error("CircularFootprintGridding is not supported on a cubed-sphere grid yet; use SubpixelGridding (quad footprints).")
    end
    metadata = InputMetadata(config)
    output_options = _output_options(config.options)
    rolling = is_rolling(time_spec)
    if rolling && backend !== nothing && compute_std
        error("Rolling-mean standard deviation is not supported on the KA/GPU backend yet; use --backend sequential, or drop --compSTD for a mean-only KA rolling mean.")
    end

    effective_n_oversample = if footprint_method isa Union{SubpixelGridding,CircularFootprintGridding} &&
                                footprint_method.n_oversample !== nothing
        footprint_method.n_oversample
    else
        n_oversample
    end

    # Snapshot schedule: spaced by `sample_step` (rolling) or `time_step` (tiled).
    snapshots = rolling ?
        collect(time_spec.start_date:time_spec.sample_step:time_spec.stop_date) :
        collect(time_spec.start_date:time_spec.time_step:time_spec.stop_date)
    n_times = length(snapshots)
    acc1, acc2 = accumulator_shape(grid_spec)
    n_vars = length(config.grid_vars)

    use_ka = backend !== nothing

    println("Output accumulator dimension (time/axis1/axis2): $n_times/$acc1/$acc2")
    rolling && println("Rolling mean: step $(time_spec.sample_step), window ±$(time_spec.window_halfwidth_days) day(s)")
    use_ka && println("Using KA backend: $backend")

    # NetCDF output is optional: an empty `outfile` means Zarr-only, which skips the
    # HDF5 write path entirely (no NetCDF byproduct, and no HDF5 failure surface).
    write_nc = !isempty(outfile)
    write_nc || zarr_out !== nothing ||
        error("Nothing to write: pass a non-empty `outfile`, `zarr_out`, or both.")

    if zarr_out !== nothing && !(grid_spec isa RectangularGridSpec)
        error("Zarr daily emit currently supports only rectangular lat/lon grids (--gridType rect); got $(typeof(grid_spec)).")
    end
    time_long_name = rolling ? "Time (UTC), center of averaging window" :
                               "Time (UTC), start of interval"
    ds_out, nc_vars, zg = nothing, nothing, nothing

    total_files = 0
    successful_files = 0
    failed_files = 0

    try
        if write_nc
            ds_out, nc_vars = create_output_dataset(outfile, grid_spec, n_times,
                config.grid_vars, compute_std; window_bounds=true,
                time_long_name=time_long_name, command=command, filters=config.filters,
                output_options=output_options)
        end
        if zarr_out !== nothing
            zg = init_zarr_store(zarr_out, grid_spec, config.grid_vars;
                                 epoch=zarr_epoch, compute_std=compute_std)
        end
        if rolling
            total_files, successful_files, failed_files = _grid_l2_rolling!(
                ds_out, nc_vars, config, grid_spec, time_spec, snapshots,
                acc1, acc2, n_vars, use_ka, backend, compute_std,
                effective_n_oversample, footprint_method, keep_going;
                zg=zg, zarr_epoch=zarr_epoch, metadata=metadata)
        else
            # Allocate work arrays (on backend for GPU, regular Array for CPU)
            if use_ka
                grid_sum = KernelAbstractions.zeros(backend, T, acc1, acc2, n_vars)
                grid_weights = KernelAbstractions.zeros(backend, T, acc1, acc2)
            else
                # Welford running mean for sequential path
                grid_data = zeros(T, acc1, acc2, n_vars)
                grid_std_arr = compute_std ? zeros(T, acc1, acc2, n_vars) : nothing
                grid_weights = zeros(T, acc1, acc2)
            end

            fill_attrib = true
            p_time = Progress(n_times, desc="Time steps: ")

            for (ct, d) in enumerate(snapshots)
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
                                                 nc_vars, fill_attrib, p_files; metadata=metadata)
                        else
                            _process_l2_file!(filepath, config, grid_spec, grid_data, grid_std_arr,
                                              grid_weights, compute_std,
                                              effective_n_oversample,
                                              footprint_method,
                                              nc_vars, fill_attrib, p_files; metadata=metadata)
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
                    _emit_l2_snapshot!(ds_out, nc_vars, zg, config, grid_spec,
                        grid_data_final, nothing, grid_w_cpu, false, ct, d;
                        window_start=d, window_end=end_date, zarr_epoch=zarr_epoch)
                    fill!(grid_sum, zero(T))
                else
                    if compute_std
                        finalize_std!(grid_std_arr, grid_weights)
                    end
                    _emit_l2_snapshot!(ds_out, nc_vars, zg, config, grid_spec,
                        grid_data, grid_std_arr, grid_weights, compute_std, ct, d;
                        window_start=d, window_end=end_date, zarr_epoch=zarr_epoch)
                    fill!(grid_data, zero(T))
                    compute_std && fill!(grid_std_arr, zero(T))
                end

                # Reset accumulators
                fill!(grid_weights, zero(T))
            end
        end

        total_files > 0 || error("No input files matched pattern '$(config.file_pattern)' in folder '$(config.folder)'")
        successful_files > 0 || error("No input files were successfully processed; failed files: $failed_files")
    finally
        _close_output_dataset!(ds_out, metadata)
    end

    if zarr_out !== nothing
        _write_zarr_provenance!(zarr_out, config, metadata, command)
        consolidate_zarr!(zarr_out)
        zarr_meta === nothing ||
            write_zarr_meta_json(zarr_meta, grid_spec, config.grid_vars; epoch=zarr_epoch)
        println("Zarr store updated: $zarr_out")
    end

    write_nc && println("Output written to: $outfile")
    nothing
end
