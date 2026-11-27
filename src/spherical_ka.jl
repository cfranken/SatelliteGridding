"""
    _footprint_cell_counts(locator, corners, normal, n, after, Val(rescan))

Count the subpoints in up to `_KA_CELL_BUFFER_SIZE` distinct cells
above `after`. The fast scan stops on overflow. The rescan retains the lowest
cell IDs, evicting the highest if necessary. An evicted ID stays above the
retention cutoff, so the next scan can count all its hits from scratch.
"""
@inline function _footprint_cell_counts(locator, corners, normal, n, after,
                                        ::Val{Rescan}) where {Rescan}
    ids = counts = _empty_ka_histogram()
    occupied = 0
    overflow = false
    @inbounds for j in 1:n
        left, right = _sphere_subpoint_edges(corners, j, n)
        for i in 1:n
            x, y, z = _sphere_subpoint(left, right, normal, i, n)
            cell = Int(locate(locator, x, y, z))
            cell <= after && continue

            ids, counts, occupied, full = _count_ka_cell(ids, counts, occupied, cell, Val(Rescan))
            if full
                overflow = true
                Rescan || return ids, counts, occupied, 0, overflow
            end
        end
    end
    return ids, counts, occupied, max(ids...), overflow
end

@inline function _spherical_footprint_geometry(vertices, fp, locator, spacing, n_override)
    @inbounds a, b, c, d = (vertices[fp, 1], vertices[fp, 2], vertices[fp, 3], vertices[fp, 4])
    corners, normal = _sphere_ordered_corners(a, b, c, d)
    n = n_override === nothing ? _unit_vector_samples(corners, spacing) : n_override
    return locator, corners, normal, n
end

# The fast scan emits nothing if the histogram overflows. Each rescan
# then counts a disjoint group of cells.
@inline function _buffer_spherical_footprint!(buffer, fp, geometry)
    nsq = geometry[4]^2
    ids, hits, occupied, _, overflow = _footprint_cell_counts(geometry..., 0, Val(false))
    after = 0
    if overflow
        ids, hits, occupied, after, overflow = _footprint_cell_counts(geometry..., after, Val(true))
    end
    while true
        @inbounds for k in 1:occupied
            push!(buffer, (fp, ids[k], hits[k], nsq))
        end
        overflow || break
        ids, hits, occupied, after, overflow = _footprint_cell_counts(geometry..., after, Val(true))
    end
    nothing
end

@kernel function _accumulate_spherical_ka!(sums, weights, pass, @Const(vertices), @Const(values),
                                           @Const(reference), locator, spacing, n_override)
    fp = @index(Global, Linear)
    geometry = _spherical_footprint_geometry(vertices, fp, locator, spacing, n_override)
    _scatter_ka_footprint!(sums, weights, pass, values, reference, fp, geometry, geometry[4]^2, Val(true))
end

function _scatter_spherical!(::CPU, sums, weights, pass, vertices, values, reference,
                             locator, spacing, n_override)
    n_fp, n_vars = size(values)
    if Threads.nthreads() == 1 || n_fp <= _KA_WORKGROUP_SIZE
        for fp in 1:n_fp
            geometry = _spherical_footprint_geometry(vertices, fp, locator, spacing, n_override)
            _scatter_ka_footprint!(sums, weights, pass, values, reference, fp, geometry,
                                   geometry[4]^2, Val(false))
        end
        return nothing
    end

    # Buffer the counted cells from contiguous ranges so the serial scatter
    # sees footprints in input order. Reuse the buffers across ranges.
    chunk_size = min(1024, cld(n_fp, Threads.nthreads()))
    n_chunks = min(Threads.nthreads(), cld(n_fp, chunk_size))
    buffers = [NTuple{4,Int}[] for _ in 1:n_chunks]
    A = eltype(weights)
    for batch in Iterators.partition(1:n_fp, chunk_size * length(buffers))
        Threads.@threads for worker in eachindex(buffers)
            buffer = empty!(buffers[worker])
            first_fp = first(batch) + (worker - 1) * chunk_size
            for fp in first_fp:min(first_fp + chunk_size - 1, last(batch))
                geometry = _spherical_footprint_geometry(vertices, fp, locator, spacing, n_override)
                _buffer_spherical_footprint!(buffer, fp, geometry)
            end
        end
        for buffer in buffers, (fp, cell, hits, nsq) in buffer
            _scatter_ka_cell!(sums, weights, pass, values, reference, fp, cell,
                              A(hits) / A(nsq), n_vars, Val(false))
        end
    end
    nothing
end

function _scatter_spherical!(backend, sums, weights, pass, vertices, values, reference,
                             locator, spacing, n_override)
    _accumulate_spherical_ka!(backend, _KA_WORKGROUP_SIZE)(
        sums, weights, pass, vertices, values, reference, locator, spacing, n_override;
        ndrange=size(values, 1))
    KernelAbstractions.synchronize(backend)
    nothing
end

"""Reusable host buffer for Cartesian corners."""
mutable struct _SphericalKAHostBatch{F}
    vertices::Matrix{NTuple{3,F}}
end

_SphericalKAHostBatch(::Type{F}) where {F} =
    _SphericalKAHostBatch(Matrix{NTuple{3,F}}(undef, 0, 0))

"""
    _accumulate_spherical_batch!(backend, sums, weights, pass, locator, F,
                                 lat, lon, values, n_override, reference, host_batch)

Upload validated corners and values and scatter one contribution per
footprint and cell. Convert degrees to Cartesian corners on the host so
cardinal angles use the same exact degree trigonometry on every backend.
The first nonempty batch sets the reference to the median of each variable's
finite values in accumulator precision, or zero if there are none. The first
pass accumulates shifted sums and weights, with bounds for runs that replay the
data. With `pass=(mean, residual)`, accumulate linear deviations into
`residual` and squared deviations into `sums` when present. Both use the cell's
first-pass mean.
"""
function _accumulate_spherical_batch!(backend, sums, weights, pass, locator, ::Type{F},
                                      lat, lon, values, n_override, reference,
                                      host_batch::_SphericalKAHostBatch{F}) where {F}
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
    device_values = _ka_batch_values!(backend, reference, values)
    spacing = _device_cell_spacing(locator, F)
    _scatter_spherical!(backend, sums, weights, pass, device_vertices, device_values,
                        reference.device, locator, spacing, n_override)
end
