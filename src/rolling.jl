# Centered temporal windows and reusable daily moments.
"""
    _grid_l2_rolling!(ds_out, nc_vars, config, grid_spec, output_options, time_spec, snapshots,
                      acc1, acc2, n_vars, use_ka, backend, compute_std,
                      effective_n_oversample, footprint_method, keep_going)
        -> (total_files, successful_files, failed_files)

Centered rolling-mean gridding. Each snapshot `c` averages the window `[c − N, c + N]`
(N = `window_halfwidth_days`). Days are gridded once into additive moments
and kept in an in-memory cache. Spherical KA merges centered daily moments
using corrected, reference-shifted means for standard deviation and sums
shifted moments for mean-only runs.
Because snapshot centers
advance monotonically, a day leaving the trailing edge is evicted and never needed again,
so each calendar day is gridded exactly once regardless of window overlap.
"""
function _grid_l2_rolling!(ds_out, nc_vars, config, grid_spec::AbstractGridSpec{T},
                           output_options, time_spec::TimeSpec, snapshots,
                           acc1, acc2, n_vars, use_ka, backend, compute_std,
                           effective_n_oversample, footprint_method, keep_going;
                           zg=nothing, zarr_epoch::Date=ZARR_EPOCH, metadata=nothing,
                           progress::Bool=true, spherical_ka=nothing) where {T}
    N = Dates.Day(time_spec.window_halfwidth_days)
    A = spherical_ka === nothing ? T : spherical_ka.accumulator_precision
    spherical_std = spherical_ka !== nothing && compute_std

    # Cache daily moments and any mean correction. Cache days without files as
    # `nothing` so we don't re-glob them across overlapping windows.
    cache = Dict{Date,Any}()
    reference = spherical_ka === nothing ? nothing : _SphericalKAReference(backend, A, n_vars)
    fill_attrib = Ref(true)
    total_files = 0
    successful_files = 0
    failed_files = 0

    # Allocate window workspaces once; finalized slices are consumed before reuse.
    win_S0 = zeros(A, acc1, acc2)
    win_S1 = zeros(A, acc1, acc2, n_vars)
    win_S2 = compute_std ? zeros(A, acc1, acc2, n_vars) : nothing
    win_correction = spherical_std ? zeros(A, acc1, acc2, n_vars) : nothing

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
                                              progress, spherical_ka, reference)
                total_files += tf
                successful_files += sf
                failed_files += ff
                cache[day] = g
            end
            dd += Dates.Day(1)
        end

        # Combine the window's days into summed moments.
        fill!(win_S0, zero(A))
        fill!(win_S1, zero(A))
        compute_std && fill!(win_S2, zero(A))
        spherical_std && fill!(win_correction, zero(A))
        dd = ws
        while dd <= we
            g = cache[Date(dd)]
            if g !== nothing
                if spherical_std
                    _merge_spherical_moments!(win_S0, win_S1, win_S2, win_correction,
                                              g.S0, g.S1, g.S2, g.correction)
                else
                    win_S0 .+= g.S0
                    win_S1 .+= g.S1
                    compute_std && (win_S2 .+= g.S2)
                end
            end
            dd += Dates.Day(1)
        end

        if spherical_ka === nothing
            finalize_moments!(win_S1, win_S2, win_S0, compute_std)
        elseif spherical_std
            _finalize_spherical_std!(win_S2, win_S0)
            _restore_spherical_reference!(win_S1, win_S0, reference, win_correction)
        else
            finalize_moments!(win_S1, nothing, win_S0, false; min_weight=zero(A))
            _restore_spherical_reference!(win_S1, win_S0, reference)
        end
        _emit_l2_snapshot!(ds_out, nc_vars, zg, config, grid_spec, output_options,
                            win_S1, win_S2, win_S0, compute_std, ct, c;
                            window_start=ws, window_end=we, zarr_epoch=zarr_epoch)

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
                  nc_vars, fill_attrib, keep_going; metadata=nothing, progress=true,
                  spherical_ka=nothing, reference=nothing)
        -> (moments, total_files, successful_files, failed_files)

Grid all input files for the single calendar day `d` into additive moments. Returns
`nothing` (no grid) when the day has no matching files. On the sequential path the
per-file Welford accumulators are converted to raw moments.
Mean-only spherical KA caches shifted sums around a shared reference.
With standard deviation, spherical KA caches the cell mean relative to that
reference as `S1` plus a low-order correction, and central M2 as `S2`.
"""
function _grid_one_day(d::DateTime, config, grid_spec::AbstractGridSpec{T},
                       acc1, acc2, n_vars, use_ka, backend, compute_std,
                       effective_n_oversample, footprint_method,
                       nc_vars, fill_attrib::Ref{Bool}, keep_going; metadata=nothing,
                       progress::Bool=true, spherical_ka=nothing, reference=nothing) where {T}
    files = find_files_range(config.file_pattern, config.folder, d, d)
    isempty(files) && return (nothing, 0, 0, 0)

    total = length(files)
    successful = 0
    failed = 0
    p_files = _progress_bar(total, "  $(Date(d)): ", progress)

    if use_ka
        A = spherical_ka === nothing ? T : spherical_ka.accumulator_precision
        day_sum = KernelAbstractions.zeros(backend, A, acc1, acc2, n_vars)
        day_w = KernelAbstractions.zeros(backend, A, acc1, acc2)
        day_m2 = compute_std ? KernelAbstractions.allocate(backend, A, size(day_sum)) : nothing
        day_residual = compute_std ? KernelAbstractions.allocate(backend, A, size(day_sum)) : nothing
        bounds = compute_std ? _spherical_mean_bounds!(day_m2, day_residual) : nothing
    else
        day_mean = zeros(T, acc1, acc2, n_vars)   # Welford running mean
        day_M2 = compute_std ? zeros(T, acc1, acc2, n_vars) : nothing     # Welford M2
        day_w = zeros(T, acc1, acc2)
    end

    std_files = String[]
    for filepath in files
        try
            if use_ka
                _process_l2_file_ka!(filepath, config, grid_spec, backend,
                                     day_sum, day_w, effective_n_oversample,
                                     footprint_method, nc_vars, fill_attrib[], p_files;
                                     metadata, spherical_ka, reference, pass=bounds)
                compute_std && push!(std_files, filepath)
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
        if compute_std
            _accumulate_l2_std_pass!(std_files, config, grid_spec, backend, day_m2, day_residual,
                day_w, day_sum, effective_n_oversample, footprint_method, spherical_ka;
                reference, bounds, progress, keep_shifted=true)
            return ((S0=Array(day_w), S1=Array(day_sum), S2=Array(day_m2),
                     correction=Array(day_residual)), total, successful, failed)
        end
        return ((S0=Array(day_w), S1=Array(day_sum), S2=nothing), total, successful, failed)
    else
        S1 = Array{T}(undef, acc1, acc2, n_vars)
        S2 = compute_std ? Array{T}(undef, acc1, acc2, n_vars) : nothing
        welford_to_moments!(S1, S2, day_mean, day_M2, day_w, compute_std)
        return ((S0=day_w, S1=S1, S2=S2), total, successful, failed)
    end
end


# Combine weighted populations using Chan's merge on the two-part means
# relative to the same run reference.
function _merge_spherical_moments!(weights, mean, m2, correction,
                                   day_weights, day_mean, day_m2, day_correction)
    @inbounds for j in axes(weights, 2), i in axes(weights, 1)
        wb = day_weights[i, j]
        wb > zero(wb) || continue
        wa = weights[i, j]
        w = wa + wb
        weights[i, j] = w
        if iszero(wa)
            for v in axes(mean, 3)
                mean[i, j, v] = day_mean[i, j, v]
                m2[i, j, v] = day_m2[i, j, v]
                correction[i, j, v] = day_correction[i, j, v]
            end
        else
            fraction = wb / w
            cross_weight = wa * fraction
            for v in axes(mean, 3)
                old_mean = mean[i, j, v]
                old_correction = correction[i, j, v]
                delta = (day_mean[i, j, v] - old_mean) + (day_correction[i, j, v] - old_correction)
                mean[i, j, v], correction[i, j, v] =
                    _spherical_two_sum(old_mean, old_correction + fraction * delta)
                m2[i, j, v] += day_m2[i, j, v] + (cross_weight * delta) * delta
            end
        end
    end
    return nothing
end
