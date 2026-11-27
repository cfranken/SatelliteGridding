"""
    grid_l2(config, grid_spec, time_spec; kwargs...)

Grid Level-2 satellite data with footprint-aware oversampling.

Reads input NetCDF files as specified by `config`, applies quality filters, subdivides
each satellite footprint into sub-pixels for proper overlap handling, and accumulates
gridded averages. Output is written to a NetCDF4 file. Quadrilateral footprints
are used by default; pass `footprint_method=CircularFootprintGridding(...)` for
circular products such as GOSAT. Circular footprints can be read from
`lat_bnd`/`lon_bnd` or from center `lat`/`lon` plus `radius`.

To archive Float64 results exactly, use Float64 accumulators and set `[output]`
`precision = "float64"`, `significant_digits_mean = "none"`, and
`significant_digits_std = "none"`. For tiled sequential runs, use Float64
bounds and spacings for `GridSpec`, or `T=Float64` for `CubedSphereGridSpec`
and `HEALPixGridSpec`. KA CPU and CUDA accumulate in Float64 by default, as do
sequential rolling runs. The CLI creates Float32 grids. A `CubedSphereGridSpec`
with `T=Float32` computes cell centre and corner coordinates in Float32 before
widening them to Float64 for NetCDF. Set `T=Float64` to compute and write those
coordinates in double precision. Use this configuration for aggregation inputs
such as HEALPix 4:1 coarsening or regridding. At fine `nside`, quantization of
nearly equal child means can dominate a small parent standard deviation.

# Arguments
- `config::DataSourceConfig`: Data source configuration (from TOML/JSON)
- `grid_spec::AbstractGridSpec`: Output grid geometry
- `time_spec::TimeSpec`: Temporal binning parameters

# Keyword Arguments
- `n_oversample=nothing`: Sub-pixel factor in `1:1024`. `nothing` uses the grid's default.
- `footprint_method::AbstractGriddingMethod=SubpixelGridding()`: Footprint geometry method.
- `compute_std::Bool=false`: Also compute per-cell standard deviation on
  sequential and KA paths. KA replays input files to refine the mean for
  Float32 accumulators and whenever standard deviation is requested.
  Mean-only Float64 KA runs use a single pass.
- `verbose::Bool=true`: Log run status with `@info`.
- `outfile::String="gridded_output.nc"`: Output file path
- `backend`: KernelAbstractions backend (default: `nothing` for sequential Welford path).
  Use `CPU()`, `CUDABackend()`, `MetalBackend()`, or another compatible KA backend.
- `geometry_precision=nothing`: Spherical KA subdivision and locator precision.
  Use `Float32` or `Float64`. Defaults to `Float64` on CUDA and KA CPU and
  `Float32` on Metal.
- `accumulator_precision=nothing`: KA precision for weights, sums, and daily moments.
  All rolling merges use Float64. Defaults to `Float64` on
  CUDA and KA CPU and `Float32` on Metal. With `Float32`, a standard deviation
  below roughly one ulp of the cell mean (`σ/|mean| ≲ 1e-7`) can be reported
  as zero. Use `Float64` when such small relative variability matters.
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
                 verbose::Bool=true,
                 command::Union{Nothing,String}=nothing)
    method = _resolve_footprint_method(footprint_method, _validate_n_oversample(n_oversample))
    run = L2Run(method, compute_std, outfile, zarr_out, zarr_meta, zarr_epoch,
                keep_going, progress, verbose, command)
    _start_l2(backend, config, grid_spec, time_spec, run, geometry_precision, accumulator_precision)
end

const _FootprintGridding = Union{SubpixelGridding,CircularFootprintGridding}

_resolve_footprint_method(method::M, n) where {M<:_FootprintGridding} =
    method.n_oversample === nothing ? M(n) : method
_resolve_footprint_method(method::AbstractGriddingMethod, n) =
    error("grid_l2 supports SubpixelGridding and CircularFootprintGridding; use grid_center for center-coordinate gridding")

"""Resolved footprint method, output destinations, and run controls for L2."""
struct L2Run{M}
    method::M
    compute_std::Bool
    outfile::String
    zarr_out::Union{Nothing,String}
    zarr_meta::Union{Nothing,String}
    zarr_epoch::Date
    keep_going::Bool
    progress::Bool
    verbose::Bool
    command::Union{Nothing,String}
end

"""Input and output state shared for every file and snapshot in an L2 run."""
struct L2Context{G,R,O,D,N,Z}
    config::DataSourceConfig
    spec::G
    run::R
    output_options::O
    ds::D
    nc_vars::N
    zg::Z
    metadata::InputMetadata
    fill_attrib::Base.RefValue{Bool}
end

_start_l2(backend::Union{Nothing,CPU}, args...) = _grid_l2(backend, args...)

# An optional backend can be loaded by resolve_backend in the calling expression.
_start_l2(backend, args...) = Base.invokelatest(_grid_l2, backend, args...)

function _l2_schedule(time_spec)
    step = time_spec.sample_step
    if step === nothing
        spacing = time_spec.time_step
        time_long_name = "Time (UTC), start of interval"
        status = nothing
    else
        spacing = step
        time_long_name = "Time (UTC), center of averaging window"
        status = "Rolling mean: step $step, window ±$(time_spec.window_halfwidth_days) day(s)"
    end
    snapshots = collect(time_spec.start_date:spacing:time_spec.stop_date)
    (; snapshots, step, time_long_name, status)
end

function _grid_l2(backend, config, spec, time_spec, run, geometry_precision, accumulator_precision)
    ka_options = _ka_options(spec, backend, geometry_precision, accumulator_precision)
    _validate(spec, run.method)
    output_options = _output_options(config.options)
    run.zarr_out === nothing || _validate_zarr_output(spec, output_options)
    _validate_l2_geometry_config(config, run.method)
    metadata = InputMetadata(config)

    (; snapshots, step, time_long_name, status) = _l2_schedule(time_spec)
    n_times = length(snapshots)
    acc1, acc2 = accumulator_shape(spec)

    if run.verbose
        @info "Output accumulator dimension (time/axis1/axis2): $n_times/$acc1/$acc2"
        status === nothing || @info status
        backend === nothing || @info "Using KA backend: $backend"
    end

    # NetCDF output is optional: an empty `outfile` means Zarr-only, which skips the
    # HDF5 write path entirely (no NetCDF byproduct, and no HDF5 failure surface).
    write_nc = !isempty(run.outfile)
    write_nc || run.zarr_out !== nothing ||
        error("Nothing to write: pass a non-empty `outfile`, `zarr_out`, or both.")

    ds, nc_vars, zg = nothing, nothing, nothing
    try
        if write_nc
            ds, nc_vars = create_output_dataset(run.outfile, spec, n_times, config.grid_vars,
                run.compute_std; window_bounds=true, time_long_name, command=run.command,
                filters=config.filters, output_options)
        end
        if run.zarr_out !== nothing
            zg = init_zarr_store(run.zarr_out, spec, config.grid_vars;
                                 epoch=run.zarr_epoch, compute_std=run.compute_std)
        end
        acc = new_accumulator(spec, backend, ka_options, length(config.grid_vars), run.compute_std, step)
        ctx = L2Context(config, spec, run, output_options, ds, nc_vars, zg, metadata, Ref(true))
        total, successful = _run_l2_schedule!(acc, ctx, time_spec, snapshots, step)
        total > 0 || error("No input files matched pattern '$(config.file_pattern)' in folder '$(config.folder)'")
        successful > 0 || error("No input files were successfully processed; failed files: $(total - successful)")
    finally
        _close_output_dataset!(ds, metadata)
    end

    if run.zarr_out !== nothing
        _write_zarr_provenance!(run.zarr_out, config, metadata, run.command)
        consolidate_zarr!(run.zarr_out)
        run.zarr_meta === nothing ||
            write_zarr_meta_json(run.zarr_meta, spec, config.grid_vars; epoch=run.zarr_epoch)
        run.verbose && @info "Zarr store updated: $(run.zarr_out)"
    end

    run.verbose && write_nc && @info "Output written to: $(run.outfile)"
    nothing
end

function _run_l2_schedule!(acc, ctx, time_spec, snapshots, ::Nothing)
    (; config, run) = ctx
    total = successful = 0
    p_time = _progress_bar(length(snapshots), "Time steps: ", run.progress)
    for (ct, d) in enumerate(snapshots)
        ProgressMeter.next!(p_time; showvalues=[(:Time, d)])
        end_date = d + time_spec.time_step * time_spec.oversample_temporal - Dates.Day(1)
        files = find_files_range(config.file_pattern, config.folder, d, end_date)
        total += length(files)

        p_files = _progress_bar(length(files), "  Files: ", run.progress)
        _reset_reference!(acc)
        successful += _accumulate_l2_files!(acc, ctx, files, p_files; keep_mean_parts=false)
        _emit_l2_snapshot!(ctx, finish!(acc), ct, d, d, end_date)
        reset!(acc)
    end
    total, successful
end

function _accumulate_l2_files!(acc, ctx, files, p_files; keep_mean_parts)
    bounds = _start_pass!(acc)
    replay_files = bounds === nothing ? nothing : String[]
    successful = 0
    for filepath in files
        try
            accumulate_file!(acc, ctx, filepath, p_files, bounds)
            replay_files === nothing || push!(replay_files, filepath)
            # A readable file with no selected observations still counts.
            successful += 1
            ctx.fill_attrib[] = false
        catch e
            if ctx.run.keep_going
                @warn "Error processing file; continuing" filepath error=sprint(showerror, e)
            else
                error("Error processing file $filepath: $(sprint(showerror, e))")
            end
        end
    end
    finish_pass!(acc, ctx, replay_files, bounds, keep_mean_parts)
    successful
end

function _emit_l2_snapshot!(ctx, arrays, ct, date, window_start, window_end)
    data, std, weights = arrays
    _emit_l2_snapshot!(ctx.ds, ctx.nc_vars, ctx.zg, ctx.config, ctx.spec, ctx.output_options,
        data, std, weights, ctx.run.compute_std, ct, date; window_start, window_end, zarr_epoch=ctx.run.zarr_epoch)
end
