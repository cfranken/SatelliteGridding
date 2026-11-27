# Centered temporal windows and reusable daily moments.

"""Daily or window weight, centered M2, and mean with optional low-order correction."""
struct DayMoments{T,M,C,S}
    weights::Matrix{T}
    mean::Array{T,3}
    m2::M
    correction::C
    overlap::S
end

"""Mean-only KA windows can add shifted sums using the common run reference."""
struct ShiftedSums{T,S}
    weights::Matrix{T}
    sums::Array{T,3}
    overlap::S
end

_host_copy(a) = Array(a)
_host_copy(::Nothing) = nothing
_window_buffer(a) = Array{Float64}(undef, size(a))
_window_buffer(::Nothing) = nothing
_window_overlap_buffer(a::AbstractMatrix{Int32}) = Matrix{Int32}(undef, size(a))
_window_overlap_buffer(a) = _window_buffer(a)

_day_moments_type(acc::WelfordAccumulator{T,M}) where {T,M} =
    DayMoments{T,M,Nothing,_host_overlap_type(acc.overlap)}
_day_moments_type(acc::AbstractKAAccumulator) = _day_moments_type(acc.moments)
_day_moments_type(acc::KAMoments{B,D,W,Nothing}) where {B,D,W} =
    ShiftedSums{eltype(W),_host_overlap_type(acc.overlap)}
_day_moments_type(acc::KAMoments) =
    _ka_day_moments_type(acc.replay, eltype(acc.weights), _host_overlap_type(acc.overlap))
_ka_day_moments_type(::KAReplay{Nothing}, ::Type{T}, ::Type{S}) where {T,S} =
    DayMoments{T,Nothing,Array{T,3},S}
_ka_day_moments_type(::KAReplay, ::Type{T}, ::Type{S}) where {T,S} =
    DayMoments{T,Array{T,3},Array{T,3},S}

day_moments(acc::WelfordAccumulator) =
    DayMoments(copy(acc.weights), copy(acc.data), _host_copy(acc.m2), nothing, _map_overlap(Array, acc.overlap))

day_moments(acc::AbstractKAAccumulator) = day_moments(acc.moments)
day_moments(acc::KAMoments{B,D,W,Nothing}) where {B,D,W} =
    ShiftedSums(Array(acc.weights), Array(acc.data), _map_overlap(Array, acc.overlap))

# Replay leaves Float64 means relative to the shared reference and Float32
# means as absolute high and low parts.
day_moments(acc::KAMoments) =
    DayMoments(Array(acc.weights), Array(acc.data), _host_copy(_ka_m2(acc.replay)),
               Array(acc.replay.residual), _map_overlap(Array, acc.overlap))

_window_moments(acc::WelfordAccumulator) = _window_moments(acc, acc.m2, nothing)
_window_moments(acc::AbstractKAAccumulator) = _window_moments(acc.moments)
_window_moments(acc::KAMoments{B,D,W,Nothing}) where {B,D,W} =
    ShiftedSums(_window_buffer(acc.weights), _window_buffer(acc.data),
                _map_overlap(_window_overlap_buffer, acc.overlap))
_window_moments(acc::KAMoments) = _window_moments(acc, _ka_m2(acc.replay), acc.replay.residual)
_window_moments(acc, m2, correction) =
    DayMoments(_window_buffer(acc.weights), _window_buffer(acc.data),
               _window_buffer(m2), _window_buffer(correction),
               _map_overlap(_window_overlap_buffer, acc.overlap))

function reset!(win::DayMoments)
    _zero!(win.weights)
    _zero!(win.mean)
    _zero!(win.m2)
    _zero!(win.correction)
    _map_overlap(_zero!, win.overlap)
    nothing
end

function reset!(win::ShiftedSums)
    _zero!(win.weights)
    _zero!(win.sums)
    _map_overlap(_zero!, win.overlap)
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
    cache = Dict{Date,Union{Nothing,_day_moments_type(acc)}}()
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
                successful += _accumulate_l2_files!(acc, ctx, files, p_files; keep_mean_parts=true)
                moments = day_moments(acc)
                reset!(acc)
                moments
            end
        end

        # Merge the cached days in calendar order.
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

function _merge_day!(win::ShiftedSums, day::ShiftedSums)
    win.weights .+= day.weights
    win.sums .+= day.sums
    _merge_overlap!(win.overlap, day.overlap)
    nothing
end

function finish!(win::ShiftedSums, acc::AbstractKAAccumulator)
    @inbounds for cell in eachindex(win.weights)
        _finalize_mean_cell!(win.sums, win.weights, cell, length(win.weights),
                            size(win.sums, 3), 0.0)
    end
    _restore_ka_reference!(win.sums, win.weights, acc.moments.reference)
    win.sums, nothing, win.weights, win.overlap
end

_finalize_rolling_mean!(win, ::WelfordAccumulator) = nothing
_finalize_rolling_mean!(win, acc::AbstractKAAccumulator) =
    _finalize_rolling_ka_mean!(win, acc.moments.reference)

function _finalize_rolling_ka_mean!(win, ::_KAReference{Float32})
    win.mean .+= win.correction
    nothing
end

_finalize_rolling_ka_mean!(win, reference::_KAReference{Float64}) =
    _restore_ka_reference!(win.mean, win.weights, reference, win.correction)

function finish!(win::DayMoments, acc)
    win.m2 === nothing || _finalize_centered_std!(win.m2, win.weights)
    _finalize_rolling_mean!(win, acc)
    win.mean, win.m2, win.weights, win.overlap
end

@inline _moment_at(::Nothing, i) = 0.0
@inline _moment_at(a, i) = @inbounds Float64(a[i])
@inline _set_moment!(::Nothing, value, i) = nothing
@inline _set_moment!(a, value, i) = @inbounds a[i] = value

@inline function _merge_mean!(mean, ::Nothing, i, old_mean, shift)
    @inbounds mean[i] = old_mean + shift
end

@inline function _merge_mean!(mean, correction, i, old_mean, shift)
    @inbounds mean[i], correction[i] = _two_sum(old_mean, shift)
end

@inline _merge_m2!(::Nothing, day_m2, i, cross_weight, delta) = nothing
@inline function _merge_m2!(m2, day_m2, i, cross_weight, delta)
    @inbounds m2[i] += Float64(day_m2[i]) + (cross_weight * delta) * delta
end

# Chan's merge. Convert both parts of each mean before arithmetic.
function _merge_day!(win::DayMoments{Float64}, day::DayMoments)
    (; weights, mean, m2, correction) = win
    stride = length(weights)
    @inbounds for v in axes(mean, 3), cell in eachindex(weights, day.weights)
        wb = Float64(day.weights[cell])
        wb > 0 || continue
        wa = weights[cell]
        i = cell + (v - 1) * stride
        if iszero(wa)
            mean[i] = day.mean[i]
            _set_moment!(m2, _moment_at(day.m2, i), i)
            _set_moment!(correction, _moment_at(day.correction, i), i)
        else
            fraction = wb / (wa + wb)
            old_mean = mean[i]
            old_correction = _moment_at(correction, i)
            delta = (Float64(day.mean[i]) - old_mean) + (_moment_at(day.correction, i) - old_correction)
            _merge_mean!(mean, correction, i, old_mean, old_correction + fraction * delta)
            _merge_m2!(m2, day.m2, i, wa * fraction, delta)
        end
    end
    # Update the shared weight after all variables have used its previous value.
    @inbounds for cell in eachindex(weights, day.weights)
        wb = Float64(day.weights[cell])
        wb > 0 && (weights[cell] += wb)
    end
    _merge_overlap!(win.overlap, day.overlap)
    nothing
end
