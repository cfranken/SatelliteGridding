# Centered temporal windows and reusable daily moments.

"""
Host arrays for one daily contribution or a rolling window. The first moment
is a sum (shifted for mean-only spherical KA) or a shifted mean (when there
is a correction). The second moment is raw for Welford and central for
spherical KA.
"""
struct DayMoments{T,M,C}
    weights::Matrix{T}
    first::Array{T,3}
    second::M
    correction::C
end

_host_buffer(a) = Array{eltype(a)}(undef, size(a))
_host_buffer(::Nothing) = nothing
_host_copy(a) = Array(a)
_host_copy(::Nothing) = nothing

function day_moments(acc::WelfordAccumulator)
    first, second = similar(acc.data), _host_buffer(acc.m2)
    welford_to_moments!(first, second, acc.data, acc.m2, acc.weights, second !== nothing)
    DayMoments(copy(acc.weights), first, second, nothing)
end

day_moments(acc::RectKAAccumulator) =
    DayMoments(Array(acc.weights), Array(acc.data), nothing, nothing)

# Replay leaves the centered moments relative to the reference shared by all days.
day_moments(acc::SphericalKAAccumulator) =
    DayMoments(Array(acc.weights), Array(acc.data), _host_copy(acc.m2), _host_copy(acc.residual))

_window_moments(acc::WelfordAccumulator) = _window_moments(acc, acc.m2, nothing)
_window_moments(acc::RectKAAccumulator) = _window_moments(acc, nothing, nothing)
_window_moments(acc::SphericalKAAccumulator) = _window_moments(acc, acc.m2, acc.residual)
_window_moments(acc, second, correction) =
    DayMoments(_host_buffer(acc.weights), _host_buffer(acc.data), _host_buffer(second), _host_buffer(correction))

function reset!(win::DayMoments)
    _zero!(win.weights)
    _zero!(win.first)
    _zero!(win.second)
    _zero!(win.correction)
    nothing
end

"""
Grid each day in the centered window at most once and merge daily moments in
calendar order. Days without files are cached as `nothing`. Because snapshot
centers advance monotonically, evicted days are never needed again.
"""
function _run_l2_schedule!(acc, ctx, time_spec, snapshots, ::Dates.Day)
    (; config, run) = ctx
    N = Dates.Day(time_spec.window_halfwidth_days)
    window = _window_moments(acc)
    cache = Dict{Date,Union{Nothing,typeof(window)}}()
    total = successful = 0

    p_time = _progress_bar(length(snapshots), "Snapshots: ", run.progress)
    for (ct, c) in enumerate(snapshots)
        ProgressMeter.next!(p_time; showvalues=[(:Center, c)])
        ws, we = c - N, c + N
        days = ws:Dates.Day(1):we

        # Grid any days in the window not already cached.
        for dd in days
            day = Date(dd)
            haskey(cache, day) && continue
            files = find_files_range(config.file_pattern, config.folder, dd, dd)
            total += length(files)
            cache[day] = if !isempty(files)
                p_files = _progress_bar(length(files), "  $day: ", run.progress)
                successful += _accumulate_l2_files!(acc, ctx, files, p_files; keep_shifted=true)
                moments = day_moments(acc)
                reset!(acc)
                moments
            end
        end

        # Combine the window's days into summed or centered moments.
        reset!(window)
        for dd in days
            moments = cache[Date(dd)]
            moments === nothing || _merge_day!(window, moments)
        end
        _emit_l2_snapshot!(ctx, finish!(window, acc), ct, c, ws, we)

        # Evict days the next (later) snapshot will no longer touch.
        if ct < length(snapshots)
            next_ws = snapshots[ct+1] - N
            filter!(entry -> DateTime(first(entry)) >= next_ws, cache)
        end
    end
    total, successful
end

function _merge_day!(win::DayMoments{T,M,Nothing}, day) where {T,M}
    win.weights .+= day.weights
    win.first .+= day.first
    win.second === nothing || (win.second .+= day.second)
    nothing
end

function _merge_day!(win::DayMoments, day)
    _merge_spherical_moments!(win.weights, win.first, win.second, win.correction,
                              day.weights, day.first, day.second, day.correction)
end

function finish!(win::DayMoments, acc)
    finalize_moments!(win.first, win.second, win.weights, win.second !== nothing)
    win.first, win.second, win.weights
end

function finish!(win::DayMoments{T,Nothing,Nothing}, acc::SphericalKAAccumulator) where {T}
    finalize_moments!(win.first, nothing, win.weights, false; min_weight=zero(T))
    _restore_spherical_reference!(win.first, win.weights, acc.reference)
    win.first, nothing, win.weights
end

function finish!(win::DayMoments, acc::SphericalKAAccumulator)
    _finalize_spherical_std!(win.second, win.weights)
    _restore_spherical_reference!(win.first, win.weights, acc.reference, win.correction)
    win.first, win.second, win.weights
end

# Combine weighted populations using Chan's merge on the two-part means
# relative to the same run reference.
function _merge_spherical_moments!(weights, mean, m2, correction,
                                   day_weights, day_mean, day_m2, day_correction)
    @inbounds for v in axes(mean, 3), j in axes(weights, 2), i in axes(weights, 1)
        wb = day_weights[i, j]
        wb > zero(wb) || continue
        wa = weights[i, j]
        if iszero(wa)
            mean[i, j, v] = day_mean[i, j, v]
            m2[i, j, v] = day_m2[i, j, v]
            correction[i, j, v] = day_correction[i, j, v]
        else
            w = wa + wb
            fraction = wb / w
            cross_weight = wa * fraction
            old_mean = mean[i, j, v]
            old_correction = correction[i, j, v]
            delta = (day_mean[i, j, v] - old_mean) + (day_correction[i, j, v] - old_correction)
            mean[i, j, v], correction[i, j, v] =
                _spherical_two_sum(old_mean, old_correction + fraction * delta)
            m2[i, j, v] += day_m2[i, j, v] + (cross_weight * delta) * delta
        end
    end
    # The previous weight is common to every variable. Update it after the merge.
    @inbounds for i in eachindex(weights, day_weights)
        wb = day_weights[i]
        wb > zero(wb) && (weights[i] += wb)
    end
    nothing
end
