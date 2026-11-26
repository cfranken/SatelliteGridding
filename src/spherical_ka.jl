# Capacity of the private histogram. Satellite footprints usually fit in one
# scan. A wide footprint uses additional scans, with the same bounded storage.
const _SPHERICAL_CELL_BUFFER_SIZE = 8
const _SPHERICAL_WORKGROUP_SIZE = 64

"""
    _footprint_cell_counts(locator, corners, normal, n, after, Val(rescan))

Count the subpoints in up to `_SPHERICAL_CELL_BUFFER_SIZE` distinct cells
above `after`. The fast scan stops on overflow. The rescan retains the lowest
cell IDs, evicting the highest if necessary. An evicted ID stays above the
retention cutoff, so the next scan can count all its hits from scratch.
"""
@inline function _footprint_cell_counts(locator, corners, normal, n, after,
                                        ::Val{Rescan}) where {Rescan}
    slots = Val(_SPHERICAL_CELL_BUFFER_SIZE)
    ids = ntuple(_ -> 0, slots)
    counts = ntuple(_ -> 0, slots)
    occupied = 0
    overflow = false
    @inbounds for j in 1:n
        left, right = _sphere_subpoint_edges(corners, j, n)
        for i in 1:n
            x, y, z = _sphere_subpoint(left, right, normal, i, n)
            cell = Int(locate(locator, x, y, z))
            cell <= after && continue

            # Fixed-index tuple comparisons keep the hot scan scalarizable.
            hit = let ids = ids
                ntuple(k -> ids[k] == cell, slots)
            end
            if any(hit)
                counts = let counts = counts
                    ntuple(k -> counts[k] + hit[k], slots)
                end
                continue
            end

            if occupied == _SPHERICAL_CELL_BUFFER_SIZE
                overflow = true
                Rescan || return ids, counts, occupied, 0, overflow
                largest = max(ids...)
                cell >= largest && continue
                slot = let ids = ids
                    sum(ntuple(k -> ifelse(ids[k] == largest, k, 0), slots))
                end
            else
                occupied += 1
                slot = occupied
            end
            ids = Base.setindex(ids, cell, slot)
            counts = Base.setindex(counts, 1, slot)
        end
    end
    return ids, counts, occupied, max(ids...), overflow
end

# Map floats to unsigned keys ordered by value so the bounds use native
# integer min/max atomics, including on Metal.
_spherical_uint(::Type{Float32}) = UInt32
_spherical_uint(::Type{Float64}) = UInt64

@inline function _spherical_ordered_key(x::A) where {A}
    U = _spherical_uint(A)
    bits = reinterpret(U, x)
    sign = one(U) << (8sizeof(U) - 1)
    return xor(bits, signbit(x) ? typemax(U) : sign)
end

@inline function _spherical_key_value(::Type{A}, key::U) where {A,U}
    sign = one(U) << (8sizeof(U) - 1)
    return reinterpret(A, xor(key, iszero(key & sign) ? typemax(U) : sign))
end

@inline _record_spherical_bounds!(::Nothing, index, x) = nothing

@inline function _record_spherical_bounds!(bounds, index, x)
    key = _spherical_ordered_key(x)
    @inbounds begin
        @atomic min(bounds.lo[index], key)
        @atomic max(bounds.hi[index], key)
    end
    return nothing
end

@inline function _scatter_spherical_counts!(sums, weights, bounds, values, reference, footprint,
                                            ids, counts, occupied, nsq)
    A = eltype(weights)
    stride = length(weights)
    @inbounds for k in 1:occupied
        cell = ids[k]
        w = A(counts[k]) / A(nsq)
        @atomic weights[cell] += w
        for v in axes(values, 2)
            index = cell + (v - 1) * stride
            x = A(values[footprint, v])
            @atomic sums[index] += w * (x - reference[v])
            _record_spherical_bounds!(bounds, index, x)
        end
    end
    return nothing
end

@inline function _scatter_spherical_counts!(m2, weights, pass::Tuple, values, reference, footprint,
                                            ids, counts, occupied, nsq)
    mean, residual = pass
    A = eltype(weights)
    stride = length(weights)
    @inbounds for k in 1:occupied
        cell = ids[k]
        w = A(counts[k]) / A(nsq)
        for v in axes(values, 2)
            index = cell + (v - 1) * stride
            delta = A(values[footprint, v]) - mean[index]
            @atomic residual[index] += w * delta
            @atomic m2[index] += w * delta * delta
        end
    end
    return nothing
end

@kernel function _accumulate_spherical_ka!(sums, weights, pass, @Const(vertices), @Const(values),
                                           @Const(reference), locator, spacing, n_override)
    footprint = @index(Global, Linear)
    @inbounds a, b, c, d = (vertices[footprint, 1], vertices[footprint, 2],
                           vertices[footprint, 3], vertices[footprint, 4])
    corners, normal = _sphere_ordered_corners(a, b, c, d)
    n = n_override === nothing ? _unit_vector_samples(corners, spacing) : n_override
    ids, counts, occupied, _, overflow = _footprint_cell_counts(locator, corners, normal,
                                                               n, 0, Val(false))
    if !overflow
        _scatter_spherical_counts!(sums, weights, pass, values, reference, footprint,
                                   ids, counts, occupied, n^2)
    else
        # The fast scan never emitted updates, so each cell receives exactly
        # one atomic update even when a footprint fills several histograms.
        after = 0
        while overflow
            ids, counts, occupied, after, overflow = _footprint_cell_counts(
                locator, corners, normal, n, after, Val(true))
            _scatter_spherical_counts!(sums, weights, pass, values, reference, footprint,
                                       ids, counts, occupied, n^2)
        end
    end
end

mutable struct _SphericalKAReference{A,D}
    host::Vector{A}
    device::D
    initialized::Bool
end

function _SphericalKAReference(backend, ::Type{T}, n_vars) where {T}
    host = Vector{T}(undef, n_vars)
    device = backend isa CPU ? host : KernelAbstractions.allocate(backend, T, size(host))
    return _SphericalKAReference(host, device, false)
end

"""Reusable host buffers for Cartesian corners and any converted science values."""
mutable struct _SphericalKAHostBatch{F,A}
    vertices::Matrix{NTuple{3,F}}
    values::Matrix{A}
end

_SphericalKAHostBatch(::Type{F}, ::Type{A}) where {F,A} =
    _SphericalKAHostBatch(Matrix{NTuple{3,F}}(undef, 0, 0), Matrix{A}(undef, 0, 0))

"""
    _accumulate_spherical_batch!(backend, sums, weights, pass, locator, F,
                                 lat, lon, values, n_override, reference, host_batch)

Upload validated corners and values and scatter one contribution per
footprint and cell. Convert degrees to Cartesian corners on the host so
cardinal angles use the same exact degree trigonometry on every backend.
The first nonempty batch sets the reference to the median of each variable's
finite values in accumulator precision, or zero if there are none. The first
pass accumulates shifted sums and weights, with
optional bounds for the standard deviation path. With `pass=(mean, residual)`
accumulate squared deviations into `sums` and linear deviations into
`residual`, both around the cell's first-pass mean.
"""
function _accumulate_spherical_batch!(backend, sums, weights, pass, locator, ::Type{F},
                                      lat, lon, values, n_override, reference,
                                      host_batch::_SphericalKAHostBatch{F,A}) where {F,A}
    n_fp = size(values, 1)
    iszero(n_fp) && return nothing
    if size(host_batch.vertices) != size(lat)
        host_batch.vertices = Matrix{NTuple{3,F}}(undef, size(lat))
    end
    vertices = host_batch.vertices
    @inbounds for j in axes(vertices, 2), i in axes(vertices, 1)
        vertices[i, j] = map(F, _sphere_unit_vector(lat[i, j], lon[i, j]))
    end
    device_vertices = _backend_input(backend, vertices)
    host_values = _spherical_host_values!(host_batch, values)
    device_values = _backend_input(backend, host_values)
    if !reference.initialized
        for v in axes(host_values, 2)
            column = filter(isfinite, @view host_values[:, v])
            reference.host[v] = isempty(column) ? zero(A) : median!(column)
        end
        reference.device === reference.host || copyto!(reference.device, reference.host)
        reference.initialized = true
    end
    spacing = _device_cell_spacing(locator, F)
    _accumulate_spherical_ka!(backend, _SPHERICAL_WORKGROUP_SIZE)(
        sums, weights, pass, device_vertices, device_values, reference.device,
        locator, spacing, n_override;
        ndrange=n_fp)
    KernelAbstractions.synchronize(backend)
    return nothing
end

_spherical_host_values!(::_SphericalKAHostBatch{F,A}, values::AbstractMatrix{A}) where {F,A} = values

function _spherical_host_values!(host_batch::_SphericalKAHostBatch{F,A}, values) where {F,A}
    if size(host_batch.values) != size(values)
        host_batch.values = Matrix{A}(undef, size(values))
    end
    host_batch.values .= A.(values)
    host_batch.values
end

# Reuse the replay buffers for pass-one bounds. Before replay, the mean
# kernel consumes the bounds and the buffers are zeroed for q and r.
function _spherical_mean_bounds!(m2, residual)
    U = _spherical_uint(eltype(m2))
    lo, hi = reinterpret(U, m2), reinterpret(U, residual)
    fill!(lo, typemax(U))
    fill!(hi, zero(U))
    return (; lo, hi)
end

@inline _bound_spherical_mean(m, ::Nothing, i) = m

@inline function _bound_spherical_mean(m::A, bounds, i) where {A}
    @inbounds lo = _spherical_key_value(A, bounds.lo[i])
    @inbounds hi = _spherical_key_value(A, bounds.hi[i])
    # A finite weighted mean is inside the observed range. This also makes
    # a constant cell exact even if its shifted sum drifted during pass one.
    return isfinite(m) ? clamp(m, lo, hi) : m
end

@kernel function _finalize_spherical_mean_ka!(sums, @Const(weights), @Const(reference), bounds)
    cell = @index(Global, Linear)
    @inbounds begin
        w = weights[cell]
        if w > zero(w)
            stride = length(weights)
            for v in axes(sums, 3)
                index = cell + (v - 1) * stride
                m = reference[v] + sums[index] / w
                sums[index] = _bound_spherical_mean(m, bounds, index)
            end
        end
    end
end

"""
    _finalize_spherical_mean!(backend, sums, weights, reference, bounds=nothing)

Convert shifted sums to absolute cell means in accumulator precision.
"""
function _finalize_spherical_mean!(backend, sums, weights, reference, bounds=nothing)
    reference.initialized || return nothing
    _finalize_spherical_mean_ka!(backend, _SPHERICAL_WORKGROUP_SIZE)(
        sums, weights, reference.device, bounds; ndrange=length(weights))
    KernelAbstractions.synchronize(backend)
    return nothing
end

# Split a sum into its rounded value and roundoff, both in accumulator precision.
@inline function _spherical_two_sum(a, b)
    s = a + b
    b_part = s - a
    return s, (a - (s - b_part)) + (b - b_part)
end

@kernel function _correct_spherical_moments_ka!(mean, m2, residual, @Const(weights),
                                               @Const(reference), ::Val{Shifted}) where {Shifted}
    cell = @index(Global, Linear)
    @inbounds begin
        w = weights[cell]
        if w > zero(w)
            stride = length(weights)
            for v in axes(mean, 3)
                index = cell + (v - 1) * stride
                r = residual[index]
                shift = r / w
                q = m2[index]
                central = q - r * shift
                refined = mean[index] + shift
                valid = isfinite(w) && isfinite(refined) && isfinite(q) &&
                        isfinite(r) && isfinite(central)
                m2[index] = valid ? max(central, zero(w)) : oftype(w, NaN)
                if Shifted
                    # Retain the subtraction roundoff too so a distant reference
                    # does not erase the low-order signal recovered by replay.
                    shifted, roundoff = _spherical_two_sum(mean[index], -reference[v])
                    mean[index] = shifted
                    residual[index] = roundoff + shift
                else
                    mean[index] = refined
                end
            end
        end
    end
end

# Preserve non-finite variance as NaN until the output writer applies fill.
function _finalize_spherical_std!(m2, weights)
    @inbounds for v in axes(m2, 3), j in axes(weights, 2), i in axes(weights, 1)
        w = weights[i, j]
        if w > zero(w)
            variance = m2[i, j, v] / w
            m2[i, j, v] = (isfinite(w) && isfinite(variance)) ? sqrt(variance) : oftype(w, NaN)
        end
    end
    return nothing
end

function _restore_spherical_reference!(mean, weights, reference, correction=nothing)
    reference.initialized || return nothing
    @inbounds for v in axes(mean, 3), j in axes(weights, 2), i in axes(weights, 1)
        weights[i, j] > 0 || continue
        if correction === nothing
            mean[i, j, v] += reference.host[v]
        else
            absolute, roundoff = _spherical_two_sum(mean[i, j, v], reference.host[v])
            mean[i, j, v] = absolute + (roundoff + correction[i, j, v])
        end
    end
    return nothing
end

_spherical_output(::Type{T}, data::Array{T}) where {T} = data
_spherical_output(::Type{T}, data::Array) where {T} = T.(data)
_spherical_output(::Type, ::Nothing) = nothing
