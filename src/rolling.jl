# Centered temporal windows and reusable daily moments.
"""
    _grid_l2_rolling!(ds_out, nc_vars, config, grid_spec, time_spec, snapshots,
                      acc1, acc2, n_vars, use_ka, backend, compute_std,
                      effective_n_oversample, footprint_method, keep_going)
        -> (total_files, successful_files, failed_files)

Centered rolling-mean gridding. Each snapshot `c` averages the window `[c − N, c + N]`
(N = `window_halfwidth_days`). Days are gridded once into additive raw moments
([`welford_to_moments!`](@ref)) and kept in an in-memory cache; because snapshot centers
advance monotonically, a day leaving the trailing edge is evicted and never needed again,
so each calendar day is gridded exactly once regardless of window overlap.
"""
function _grid_l2_rolling!(ds_out, nc_vars, config, grid_spec::AbstractGridSpec{T},
                           time_spec::TimeSpec, snapshots,
                           acc1, acc2, n_vars, use_ka, backend, compute_std,
                           effective_n_oversample, footprint_method, keep_going;
                           zg=nothing, zarr_epoch::Date=ZARR_EPOCH, metadata=nothing,
                           progress::Bool=true) where {T}
    N = Dates.Day(time_spec.window_halfwidth_days)

    # day (calendar Date) => (S0, S1, S2) raw moments, or `nothing` when the day had no
    # files (cached as a negative result so we don't re-glob it across overlapping windows).
    cache = Dict{Date,Any}()
    fill_attrib = Ref(true)
    total_files = 0
    successful_files = 0
    failed_files = 0

    # Allocate window workspaces once; finalized slices are consumed before reuse.
    win_S0 = zeros(T, acc1, acc2)
    win_S1 = zeros(T, acc1, acc2, n_vars)
    win_S2 = compute_std ? zeros(T, acc1, acc2, n_vars) : nothing

    p_time = _progress_bar(length(snapshots), "Snapshots: ", progress)
    for (ct, c) in enumerate(snapshots)
        ProgressMeter.next!(p_time; showvalues=[(:Center, c)])
        ws = c - N
        we = c + N

        # Grid any days in the window not already cached.
        dd = ws
        while dd <= we
            day = Date(dd)
            if !haskey(cache, day)
                g, tf, sf, ff = _grid_one_day(dd, config, grid_spec, acc1, acc2, n_vars,
                                              use_ka, backend, compute_std,
                                              effective_n_oversample, footprint_method,
                                              nc_vars, fill_attrib, keep_going; metadata=metadata,
                                              progress=progress)
                total_files += tf
                successful_files += sf
                failed_files += ff
                cache[day] = g
            end
            dd += Dates.Day(1)
        end

        # Combine the window's days into summed moments.
        fill!(win_S0, zero(T))
        fill!(win_S1, zero(T))
        compute_std && fill!(win_S2, zero(T))
        dd = ws
        while dd <= we
            g = cache[Date(dd)]
            if g !== nothing
                win_S0 .+= g.S0
                win_S1 .+= g.S1
                compute_std && (win_S2 .+= g.S2)
            end
            dd += Dates.Day(1)
        end

        finalize_moments!(win_S1, win_S2, win_S0, compute_std)
        _emit_l2_snapshot!(ds_out, nc_vars, zg, config, grid_spec, win_S1, win_S2,
                            win_S0, compute_std, ct, c; window_start=ws, window_end=we,
                            zarr_epoch=zarr_epoch)

        # Evict days the next (later) snapshot will no longer touch.
        if ct < length(snapshots)
            next_ws = snapshots[ct+1] - N
            for day in collect(keys(cache))
                DateTime(day) < next_ws && delete!(cache, day)
            end
        end
    end

    return total_files, successful_files, failed_files
end

"""
    _grid_one_day(d, config, grid_spec, acc1, acc2, n_vars, use_ka, backend,
                  compute_std, effective_n_oversample, footprint_method,
                  nc_vars, fill_attrib, keep_going; metadata=nothing, progress=true)
        -> ((S0, S1, S2) | nothing, total_files, successful_files, failed_files)

Grid all input files for the single calendar day `d` into additive raw moments. Returns
`nothing` (no grid) when the day has no matching files. On the sequential path the
per-file Welford accumulators are converted to moments; on the KA path the accumulators
are already weighted sums (`S0 = Σw`, `S1 = Σw·x`).
"""
function _grid_one_day(d::DateTime, config, grid_spec::AbstractGridSpec{T},
                       acc1, acc2, n_vars, use_ka, backend, compute_std,
                       effective_n_oversample, footprint_method,
                       nc_vars, fill_attrib::Ref{Bool}, keep_going; metadata=nothing,
                       progress::Bool=true) where {T}
    files = find_files_range(config.file_pattern, config.folder, d, d)
    isempty(files) && return (nothing, 0, 0, 0)

    total = length(files)
    successful = 0
    failed = 0
    p_files = _progress_bar(total, "  $(Date(d)): ", progress)

    if use_ka
        day_sum = KernelAbstractions.zeros(backend, T, acc1, acc2, n_vars)
        day_w = KernelAbstractions.zeros(backend, T, acc1, acc2)
    else
        day_mean = zeros(T, acc1, acc2, n_vars)   # Welford running mean
        day_M2 = compute_std ? zeros(T, acc1, acc2, n_vars) : nothing     # Welford M2
        day_w = zeros(T, acc1, acc2)
    end

    for filepath in files
        try
            if use_ka
                _process_l2_file_ka!(filepath, config, grid_spec, backend,
                                     day_sum, day_w, effective_n_oversample,
                                     footprint_method, nc_vars, fill_attrib[], p_files; metadata=metadata)
            else
                _process_l2_file!(filepath, config, grid_spec, day_mean, day_M2,
                                  day_w, compute_std, effective_n_oversample,
                                  footprint_method, nc_vars, fill_attrib[], p_files; metadata=metadata)
            end
            successful += 1
            fill_attrib[] = false
        catch e
            failed += 1
            if keep_going
                @warn "Error processing file; continuing" filepath error=sprint(showerror, e)
            else
                error("Error processing file $filepath: $(sprint(showerror, e))")
            end
        end
    end

    if use_ka
        # KA accumulators are already weighted sums; std unsupported on this path.
        return ((S0=Array(day_w), S1=Array(day_sum), S2=nothing), total, successful, failed)
    else
        S1 = Array{T}(undef, acc1, acc2, n_vars)
        S2 = compute_std ? Array{T}(undef, acc1, acc2, n_vars) : nothing
        welford_to_moments!(S1, S2, day_mean, day_M2, day_w, compute_std)
        return ((S0=day_w, S1=S1, S2=S2), total, successful, failed)
    end
end
