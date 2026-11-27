# Bound private kernel storage. Wide footprints rescan above the largest cell
# retained in the previous scan.
const _KA_CELL_BUFFER_SIZE = 8
const _KA_WORKGROUP_SIZE = 64

@inline _empty_ka_histogram() = ntuple(_ -> 0, Val(_KA_CELL_BUFFER_SIZE))

@inline function _count_ka_cell(ids, counts, occupied, cell, ::Val{Rescan}) where {Rescan}
    slots = Val(_KA_CELL_BUFFER_SIZE)
    # Fixed-index tuple comparisons keep the hot scan scalarizable.
    hit = ntuple(k -> ids[k] == cell, slots)
    if any(hit)
        return ids, ntuple(k -> counts[k] + hit[k], slots), occupied, false
    end
    full = occupied == _KA_CELL_BUFFER_SIZE
    if full
        Rescan || return ids, counts, occupied, true
        largest = max(ids...)
        cell >= largest && return ids, counts, occupied, true
        slot = sum(ntuple(k -> ifelse(ids[k] == largest, k, 0), slots))
    else
        occupied += 1
        slot = occupied
    end
    return Base.setindex(ids, cell, slot), Base.setindex(counts, 1, slot), occupied, full
end

# Map floats to unsigned keys ordered by value so the bounds use native
# integer min/max atomics, including on Metal.
_ka_uint(::Type{Float32}) = UInt32
_ka_uint(::Type{Float64}) = UInt64

@inline function _ka_ordered_key(x::A) where {A}
    U = _ka_uint(A)
    bits = reinterpret(U, x)
    sign = one(U) << (8sizeof(U) - 1)
    return xor(bits, signbit(x) ? typemax(U) : sign)
end

@inline function _ka_key_value(::Type{A}, key::U) where {A,U}
    sign = one(U) << (8sizeof(U) - 1)
    return reinterpret(A, xor(key, iszero(key & sign) ? typemax(U) : sign))
end

@inline _record_ka_bounds!(::Nothing, index, x, atomic) = nothing

@inline function _record_ka_bounds!(bounds::NamedTuple, index, x, ::Val{false})
    key = _ka_ordered_key(x)
    @inbounds begin
        bounds.lo[index] = min(bounds.lo[index], key)
        bounds.hi[index] = max(bounds.hi[index], key)
    end
    return nothing
end

@inline function _record_ka_bounds!(bounds::NamedTuple, index, x, ::Val{true})
    key = _ka_ordered_key(x)
    @inbounds begin
        @atomic min(bounds.lo[index], key)
        @atomic max(bounds.hi[index], key)
    end
    return nothing
end

@inline function _add_ka_cell!(a, index, x, ::Val{true})
    @inbounds @atomic a[index] += x
    nothing
end

@inline function _add_ka_cell!(a, index, x, ::Val{false})
    @inbounds a[index] += x
    nothing
end

@inline function _scatter_ka_cell!(sums, weights, bounds, values, reference,
                                   footprint, cell, w, n_vars, atomic)
    A = eltype(weights)
    stride = length(weights)
    @inbounds begin
        _add_ka_cell!(weights, cell, w, atomic)
        for v in 1:n_vars
            index = cell + (v - 1) * stride
            x = A(values[footprint, v])
            _add_ka_cell!(sums, index, w * (x - reference[v]), atomic)
            _record_ka_bounds!(bounds, index, x, atomic)
        end
    end
    return nothing
end

@inline _add_ka_squared_deviation!(::Nothing, index, w, delta, atomic) = nothing

@inline _add_ka_squared_deviation!(m2, index, w, delta, atomic) =
    _add_ka_cell!(m2, index, w * delta * delta, atomic)

@inline function _scatter_ka_cell!(m2, weights, pass::Tuple, values, reference,
                                   footprint, cell, w, n_vars, atomic)
    mean, residual = pass
    A = eltype(weights)
    stride = length(weights)
    @inbounds for v in 1:n_vars
        index = cell + (v - 1) * stride
        delta = A(values[footprint, v]) - mean[index]
        _add_ka_cell!(residual, index, w * delta, atomic)
        _add_ka_squared_deviation!(m2, index, w, delta, atomic)
    end
    return nothing
end

@inline function _scatter_ka_counts!(sums, weights, pass, values, reference, footprint,
                                     ids, counts, occupied, denominator, n_vars, atomic)
    A = eltype(weights)
    @inbounds for k in 1:occupied
        w = A(counts[k]) / A(denominator)
        _scatter_ka_cell!(sums, weights, pass, values, reference,
                          footprint, ids[k], w, n_vars, atomic)
    end
    return nothing
end

# The fast scan emits nothing if the histogram overflows. Each rescan then
# counts a disjoint group of cells and sends exactly one update per cell.
@inline function _scatter_ka_footprint!(sums, weights, pass, values, reference, fp,
                                       geometry, denominator, atomic)
    ids, hits, occupied, _, overflow = _footprint_cell_counts(geometry..., 0, Val(false))
    if !overflow
        _scatter_ka_counts!(sums, weights, pass, values, reference, fp,
                            ids, hits, occupied, denominator, size(values, 2), atomic)
    else
        after = 0
        while overflow
            ids, hits, occupied, after, overflow = _footprint_cell_counts(geometry..., after, Val(true))
            _scatter_ka_counts!(sums, weights, pass, values, reference, fp,
                                ids, hits, occupied, denominator, size(values, 2), atomic)
        end
    end
    nothing
end

mutable struct _KAReference{A,D}
    host::Vector{A}
    device::D
    values::Matrix{A}
    initialized::Bool
end

function _KAReference(backend, ::Type{T}, n_vars) where {T}
    host = Vector{T}(undef, n_vars)
    device = backend isa CPU ? host : KernelAbstractions.allocate(backend, T, size(host))
    return _KAReference(host, device, Matrix{T}(undef, 0, 0), false)
end

function _ka_batch_values!(backend, reference::_KAReference{A}, values) where {A}
    host_values = _ka_host_values!(reference, values)
    device_values = _backend_input(backend, host_values)
    if !reference.initialized
        for v in axes(host_values, 2)
            column = filter(isfinite, @view host_values[:, v])
            reference.host[v] = isempty(column) ? zero(A) : median!(column)
        end
        reference.device === reference.host || copyto!(reference.device, reference.host)
        reference.initialized = true
    end
    device_values
end

_ka_host_values!(::_KAReference{A}, values::AbstractMatrix{A}) where {A} = values

function _ka_host_values!(reference::_KAReference{A}, values) where {A}
    if size(reference.values) != size(values)
        reference.values = Matrix{A}(undef, size(values))
    end
    reference.values .= A.(values)
    reference.values
end

function _reset_ka_bounds!(bounds)
    U = eltype(bounds.lo)
    fill!(bounds.lo, typemax(U))
    fill!(bounds.hi, zero(U))
    return bounds
end

@inline _bound_ka_mean(m, ::Nothing, i) = m

@inline function _ka_cell_bounds(::Type{A}, bounds, i) where {A}
    @inbounds lo = _ka_key_value(A, bounds.lo[i])
    @inbounds hi = _ka_key_value(A, bounds.hi[i])
    return lo, hi
end

@inline function _bound_ka_mean(m::A, bounds, i) where {A}
    # A finite weighted mean is inside the observed range. This also makes
    # a constant cell exact even if its shifted sum drifted during pass one.
    lo, hi = _ka_cell_bounds(A, bounds, i)
    return isfinite(m) ? clamp(m, lo, hi) : m
end

@kernel function _finalize_ka_mean_ka!(sums, @Const(weights), @Const(reference), bounds)
    cell = @index(Global, Linear)
    @inbounds begin
        w = weights[cell]
        if w > zero(w)
            stride = length(weights)
            for v in axes(sums, 3)
                index = cell + (v - 1) * stride
                m = reference[v] + sums[index] / w
                sums[index] = _bound_ka_mean(m, bounds, index)
            end
        end
    end
end

"""
    _finalize_ka_mean!(backend, sums, weights, reference, bounds=nothing)

Convert shifted sums to absolute cell means in accumulator precision.
"""
function _finalize_ka_mean!(backend, sums, weights, reference, bounds=nothing)
    reference.initialized || return nothing
    _finalize_ka_mean_ka!(backend, _KA_WORKGROUP_SIZE)(
        sums, weights, reference.device, bounds; ndrange=length(weights))
    KernelAbstractions.synchronize(backend)
    return nothing
end

# Split a sum into its rounded value and roundoff, both in accumulator precision.
@inline function _two_sum(a, b)
    s = a + b
    b_part = s - a
    return s, (a - (s - b_part)) + (b - b_part)
end

@inline _refine_ka_mean(m, shift, ::Nothing, index) = (m + shift, m, shift)

@inline function _refine_ka_mean(m::A, shift, bounds, index) where {A}
    refined, low = _two_sum(m, shift)
    if isfinite(refined)
        lo, hi = _ka_cell_bounds(A, bounds, index)
        bounded = clamp(refined, lo, hi)
        # The low part can move a rounded endpoint outside the range when
        # the two parts are promoted to Float64 for a rolling merge.
        if refined != bounded || (refined == lo && low < 0) || (refined == hi && low > 0)
            return bounded, bounded, zero(A)
        end
    end
    return refined, refined, low
end

@inline _ka_mean_parts(base, correction, reference, bounds) = (base, correction)

@inline function _ka_mean_parts(base, correction, reference, ::Nothing)
    # Float64 rolling days share a reference. Preserve the subtraction
    # roundoff so that a distant reference retains the replay correction.
    shifted, roundoff = _two_sum(base, -reference)
    return shifted, roundoff + correction
end

@inline _center_ka_m2!(::Nothing, index, w, r, shift, refined) = nothing

@inline function _center_ka_m2!(m2, index, w, r, shift, refined)
    @inbounds begin
        q = m2[index]
        central = q - r * shift
        valid = isfinite(w) && isfinite(refined) && isfinite(q) &&
                isfinite(r) && isfinite(central)
        m2[index] = valid ? max(central, zero(w)) : oftype(w, NaN)
    end
    nothing
end

@kernel function _correct_ka_moments_ka!(mean, m2, residual, @Const(weights),
                                       @Const(reference), @Const(bounds), ::Val{Rolling}) where {Rolling}
    cell = @index(Global, Linear)
    @inbounds begin
        w = weights[cell]
        if w > zero(w)
            stride = length(weights)
            for v in axes(mean, 3)
                index = cell + (v - 1) * stride
                r = residual[index]
                shift = r / w
                refined, base, correction = _refine_ka_mean(mean[index], shift, bounds, index)
                _center_ka_m2!(m2, index, w, r, shift, refined)
                if Rolling
                    mean[index], residual[index] = _ka_mean_parts(base, correction, reference[v], bounds)
                else
                    mean[index] = refined
                end
            end
        end
    end
end

# Preserve non-finite variance as NaN until the output writer applies fill.
function _finalize_centered_std!(m2, weights)
    @inbounds for v in axes(m2, 3), j in axes(weights, 2), i in axes(weights, 1)
        w = weights[i, j]
        if w > zero(w)
            variance = m2[i, j, v] / w
            m2[i, j, v] = isfinite(w) && isfinite(variance) ?
                sqrt(max(variance, zero(variance))) : oftype(w, NaN)
        end
    end
    return nothing
end

function _restore_ka_reference!(mean, weights, reference, correction=nothing)
    reference.initialized || return nothing
    @inbounds for v in axes(mean, 3), j in axes(weights, 2), i in axes(weights, 1)
        weights[i, j] > 0 || continue
        if correction === nothing
            mean[i, j, v] += reference.host[v]
        else
            absolute, roundoff = _two_sum(mean[i, j, v], reference.host[v])
            mean[i, j, v] = absolute + (roundoff + correction[i, j, v])
        end
    end
    return nothing
end
