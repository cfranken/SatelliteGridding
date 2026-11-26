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
- `grid_spec::AbstractGridSpec`: Output grid geometry
- `time_spec::TimeSpec`: Temporal binning parameters

# Keyword Arguments
- `n_oversample=nothing`: Sub-pixel factor in `1:1024`. `nothing` uses the grid's default.
- `footprint_method::AbstractGriddingMethod=SubpixelGridding()`: Footprint geometry method.
- `compute_std::Bool=false`: Also compute per-cell standard deviation. Mean-only
  spherical KA runs skip mean refinement. For accumulator precision `A`, values
  with magnitude below about `eps(A) * abs(reference - value)` can round away
  when shifted by the per-variable reference (`eps(Float64)` ≈ `2e-16`,
  `eps(Float32)` ≈ `1e-7`). Enable `compute_std` to recover those values.
- `outfile::String="gridded_output.nc"`: Output file path
- `backend`: KernelAbstractions backend (default: `nothing` for sequential Welford path).
  Use `CPU()`, `CUDABackend()`, `MetalBackend()`, or another compatible KA backend.
- `geometry_precision=nothing`: Spherical KA subdivision and locator precision.
  Use `Float32` or `Float64`. Defaults to `Float64` on CUDA and KA CPU and
  `Float32` on Metal.
- `accumulator_precision=nothing`: Spherical KA precision for weights, sums,
  and variance, including rolling daily moments. Defaults to `Float64` on
  CUDA and KA CPU and `Float32` on Metal. With `Float32`, a standard deviation
  below roughly one ulp of the cell mean (`σ/|mean| ≲ 1e-7`) can be reported
  as zero. Use `Float64` when such small relative variability matters.
  The grid element type is used for finalized slices.
- `keep_going::Bool=false`: Continue after per-file errors. By default, the
  first file error stops the run.
- `progress::Bool=stderr isa Base.TTY`: Draw progress bars (on stderr). By default they're
  drawn only when stderr is a terminal: redirected to a file or a pipe, as in batch jobs, every
  redraw would become new lines.
"""
function grid_l2(config::DataSourceConfig, grid_spec::AbstractGridSpec, time_spec::TimeSpec;
                 n_oversample=nothing,
                 footprint_method::AbstractGriddingMethod=SubpixelGridding(),
                 compute_std::Bool=false,
                 outfile::String="gridded_output.nc",
                 zarr_out::Union{Nothing,String}=nothing,
                 zarr_meta::Union{Nothing,String}=nothing,
                 zarr_epoch::Date=ZARR_EPOCH,
                 backend=nothing,
                 geometry_precision=nothing,
                 accumulator_precision=nothing,
                 keep_going::Bool=false,
                 progress::Bool=_default_progress(),
                 command::Union{Nothing,String}=nothing)
    n_oversample = _validate_n_oversample(n_oversample)
    footprint_method isa SubpixelGridding ||
        footprint_method isa CircularFootprintGridding ||
        error("grid_l2 supports SubpixelGridding and CircularFootprintGridding; use grid_center for center-coordinate gridding")
    method_n_oversample = footprint_method.n_oversample
    effective_n_oversample = method_n_oversample === nothing ? n_oversample : method_n_oversample
    options = (; effective_n_oversample, footprint_method, compute_std, outfile, zarr_out,
               zarr_meta, zarr_epoch, backend, geometry_precision, accumulator_precision,
               keep_going, progress, command)
    if backend === nothing || backend isa CPU
        return _grid_l2(config, grid_spec, time_spec; options...)
    end
    # An optional backend can be loaded by resolve_backend in the calling expression.
    return Base.invokelatest(_grid_l2, config, grid_spec, time_spec; options...)
end

function _grid_l2(config::DataSourceConfig, grid_spec::AbstractGridSpec{T}, time_spec::TimeSpec;
                  effective_n_oversample, footprint_method, compute_std, outfile, zarr_out, zarr_meta,
                  zarr_epoch, backend, geometry_precision, accumulator_precision,
                  keep_going, progress, command) where {T}

    spherical_ka = _spherical_ka_options(grid_spec, backend, geometry_precision, accumulator_precision)
    A = spherical_ka === nothing ? T : spherical_ka.accumulator_precision
    spherical_std = spherical_ka !== nothing && compute_std
    grid_spec isa AbstractSphericalGridSpec && footprint_method isa CircularFootprintGridding &&
        error("CircularFootprintGridding is unsupported for $(typeof(grid_spec)). Use SubpixelGridding (quad footprints).")
    zarr_out === nothing || _validate_zarr_grid(grid_spec)
    _validate_l2_geometry_config(config, footprint_method)
    metadata = InputMetadata(config)
    output_options = _output_options(config.options)
    rolling = is_rolling(time_spec)
    if rolling && backend !== nothing && compute_std && !spherical_std
        error("Rolling-mean standard deviation is not supported on the KA/GPU backend yet; use --backend sequential, or drop --compSTD for a mean-only KA rolling mean.")
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
                output_options=output_options, spherical_ka_std=spherical_std)
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
                zg=zg, zarr_epoch=zarr_epoch, metadata=metadata, progress=progress, spherical_ka)
        else
            # Allocate work arrays (on backend for GPU, regular Array for CPU)
            if use_ka
                grid_sum = KernelAbstractions.zeros(backend, A, acc1, acc2, n_vars)
                grid_weights = KernelAbstractions.zeros(backend, A, acc1, acc2)
                grid_m2 = spherical_std ? KernelAbstractions.allocate(backend, A, size(grid_sum)) : nothing
                grid_residual = spherical_std ? KernelAbstractions.allocate(backend, A, size(grid_sum)) : nothing
            else
                # Welford running mean for sequential path
                grid_data = zeros(T, acc1, acc2, n_vars)
                grid_std_arr = compute_std ? zeros(T, acc1, acc2, n_vars) : nothing
                grid_weights = zeros(T, acc1, acc2)
            end

            reference = spherical_ka === nothing ? nothing : _SphericalKAReference(backend, A, n_vars)
            fill_attrib = true
            p_time = _progress_bar(n_times, "Time steps: ", progress)

            for (ct, d) in enumerate(snapshots)
                ProgressMeter.next!(p_time; showvalues=[(:Time, d)])

                # Find files for this time window
                end_date = d + time_spec.time_step * time_spec.oversample_temporal - Dates.Day(1)
                files = find_files_range(config.file_pattern, config.folder, d, end_date)
                total_files += length(files)

                p_files = _progress_bar(length(files), "  Files: ", progress)

                reference === nothing || (reference.initialized = false)
                bounds = spherical_std ? _spherical_mean_bounds!(grid_m2, grid_residual) : nothing
                std_files = String[]
                for filepath in files
                    try
                        if use_ka
                            _process_l2_file_ka!(filepath, config, grid_spec, backend,
                                                 grid_sum, grid_weights,
                                                 effective_n_oversample,
                                                 footprint_method,
                                                 nc_vars, fill_attrib, p_files;
                                                 metadata, spherical_ka, reference, pass=bounds)
                            spherical_std && push!(std_files, filepath)
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
                    if spherical_std
                        _accumulate_l2_std_pass!(std_files, config, grid_spec, backend,
                            grid_m2, grid_residual, grid_weights, grid_sum, effective_n_oversample,
                            footprint_method, spherical_ka; reference, bounds, progress)
                    elseif spherical_ka !== nothing
                        _finalize_spherical_mean!(backend, grid_sum, grid_weights, reference)
                    end
                    grid_data_final = Array(grid_sum)
                    grid_w_cpu = Array(grid_weights)
                    std_final = if spherical_std
                        m2 = Array(grid_m2)
                        _finalize_spherical_std!(m2, grid_w_cpu)
                        m2
                    else
                        spherical_ka === nothing && finalize_mean!(grid_data_final, grid_w_cpu)
                        nothing
                    end
                    _emit_l2_snapshot!(ds_out, nc_vars, zg, config, grid_spec,
                        _spherical_output(T, grid_data_final), _spherical_output(T, std_final),
                        _spherical_output(T, grid_w_cpu), spherical_std, ct, d;
                        window_start=d, window_end=end_date, zarr_epoch=zarr_epoch)
                    fill!(grid_sum, zero(A))
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
                fill!(grid_weights, zero(A))
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
