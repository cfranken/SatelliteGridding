"""
    accumulate_footprint!(grid_data, grid_std, grid_weights, compute_std,
                          lat_idx, lon_idx, values, n_pixels, n_vars, n_oversample,
                          points_buf, ix_buf, iy_buf,
                          lats0_buf, lons0_buf, lats1_buf, lons1_buf)

Accumulate a batch of satellite observations into the output grid using footprint-aware
oversampling and Welford's online averaging algorithm.

For each pixel:
- **Fast path**: If all 4 corner coordinates fall into the same grid cell, the full
  observation value is assigned to that cell (weight=1).
- **Oversample path**: If the footprint spans multiple cells, it is subdivided into
  `n_oversample × n_oversample` sub-pixels. Each sub-pixel contributes weight `1/n²`
  to whichever grid cell it falls in.

Sub-pixels (or single-cell footprints) that fall outside `1..n_lon × 1..n_lat` are
dropped without renormalizing remaining weights so boundary-cell weights reflect true
fractional footprint overlap and remain domain-invariant. On global rectangular grids,
samples beyond `±180°` longitude or beyond the poles are dropped rather than wrapped.

The Welford algorithm maintains a running mean and sum-of-squared-deviations (M2),
enabling numerically stable computation of mean and variance in a single pass.

# Arguments
- `grid_data`: Output array `(n_lon, n_lat, n_vars)` — running mean
- `grid_std`: Output array `(n_lon, n_lat, n_vars)` — running M2 (sum of squared deviations); may be `nothing` when `compute_std=false`
- `grid_weights`: Output array `(n_lon, n_lat)` — accumulated weights
- `compute_std`: Whether to track variance
- `lat_idx`, `lon_idx`: Pre-computed fractional grid indices for corners, size `(n_pixels, 4)`
- `values`: Input data, size `(n_pixels, n_vars)`
- `n_pixels`: Number of observations in this batch
- `n_vars`: Number of variables being gridded
- `n_oversample`: Sub-pixel subdivision factor (n)
- Remaining arguments are pre-allocated work buffers
"""
function accumulate_footprint!(grid_data::AbstractArray{T,3},
                               grid_std::Union{Nothing,AbstractArray{T,3}},
                               grid_weights::AbstractMatrix{T},
                               compute_std::Bool,
                               lat_idx::AbstractMatrix, lon_idx::AbstractMatrix,
                               values::AbstractMatrix,
                               n_pixels::Int, n_vars::Int, n_oversample::Int,
                               points_buf::AbstractArray{<:AbstractFloat,3},
                               ix_buf::AbstractVector{<:Integer},
                               iy_buf::AbstractVector{<:Integer},
                               lats0_buf::AbstractVector,
                               lons0_buf::AbstractVector,
                               lats1_buf::AbstractVector,
                               lons1_buf::AbstractVector) where {T}
    n_lon, n_lat = size(grid_weights)

    # Compute grid-cell indices for all 4 corners
    ilon = floor.(Int32, lon_idx)
    ilat = floor.(Int32, lat_idx)
    min_lat = minimum(ilat, dims=2)
    max_lat = maximum(ilat, dims=2)
    min_lon = minimum(ilon, dims=2)
    max_lon = maximum(ilon, dims=2)
    dim_lat = max_lat .- min_lat
    dim_lon = max_lon .- min_lon

    fac = T(1 / n_oversample^2)

    @inbounds for i in 1:n_pixels
        if (dim_lat[i] == 0) & (dim_lon[i] == 0)
            # Fast path: all corners in one cell
            col = ilon[i, 1]
            row = ilat[i, 1]
            if 1 <= col <= n_lon && 1 <= row <= n_lat
                _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                                 col, row, values, i, n_vars, T(1))
            end
        elseif dim_lon[i] < n_oversample
            # Oversample path: footprint spans multiple cells
            compute_subpixels!(points_buf, @view(lat_idx[i, :]), @view(lon_idx[i, :]), n_oversample,
                               lats0_buf, lons0_buf, lats1_buf, lons1_buf)
            floor_indices!(ix_buf, iy_buf, points_buf)

            @inbounds for j in eachindex(ix_buf)
                col = iy_buf[j]
                row = ix_buf[j]
                if 1 <= col <= n_lon && 1 <= row <= n_lat
                    _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                                     col, row, values, i, n_vars, fac)
                end
            end
        end
    end
    nothing
end

"""
    accumulate_footprint_cs!(grid_data, grid_std, grid_weights, compute_std, spec,
                             lat_c, lon_c, values, n_pixels, n_vars, n_oversample,
                             points_buf, lats0_buf, lons0_buf, lats1_buf, lons1_buf)

Footprint-aware accumulation onto a [`CubedSphereGridSpec`](@ref) folded into the
`(Nc, 6·Nc)` accumulator. Unlike the rectangular [`accumulate_footprint!`], the
subdivision happens in **lon/lat space** because the cubed-sphere map is nonlinear
and discontinuous across panel boundaries. Each subpoint is then projected
independently via [`to_fractional_index`], so a footprint straddling a cube edge
correctly deposits weight on both panels.

For each pixel:
- **Fast path**: all 4 corners project to the same `(col, row)` cell → full value, weight 1.
- **Oversample path**: unwrap corner longitudes relative to the first corner
  before subdividing the quadrilateral into `n × n` subpoints in lon/lat space.
  Project each subpoint and add weight `1/n²` to the cell it lands in.

`lat_c` and `lon_c` are `(n_pixels, 4)` corner coordinates in degrees, sorted
counterclockwise. Buffers match [`compute_subpixels!`]. HEALPix uses a separate
unit-vector subdivision path.
"""
function accumulate_footprint_cs!(grid_data, grid_std, grid_weights, compute_std,
                                  spec::CubedSphereGridSpec, args...)
    _accumulate_spherical_footprint!(grid_data, grid_std, grid_weights, compute_std,
                                    spec, args...)
end

function _accumulate_spherical_footprint!(grid_data::AbstractArray{T,3},
                                          grid_std::Union{Nothing,AbstractArray{T,3}},
                                          grid_weights::AbstractMatrix{T},
                                          compute_std::Bool,
                                          spec::AbstractSphericalGridSpec,
                                          lat_c::AbstractMatrix, lon_c::AbstractMatrix,
                                          values::AbstractMatrix,
                                          n_pixels::Int, n_vars::Int, n_oversample::Int,
                                          points_buf::AbstractArray{<:AbstractFloat,3},
                                          lats0_buf::AbstractVector, lons0_buf::AbstractVector,
                                          lats1_buf::AbstractVector, lons1_buf::AbstractVector) where {T}
    fac = T(1 / n_oversample^2)
    lonv = Vector{Float64}(undef, 4)
    latv = Vector{Float64}(undef, 4)

    @inbounds for i in 1:n_pixels
        # Project the 4 corners; fast-path when they all share one cell.
        c1, r1 = _spherical_cell(spec, lat_c[i, 1], lon_c[i, 1])
        c2, r2 = _spherical_cell(spec, lat_c[i, 2], lon_c[i, 2])
        c3, r3 = _spherical_cell(spec, lat_c[i, 3], lon_c[i, 3])
        c4, r4 = _spherical_cell(spec, lat_c[i, 4], lon_c[i, 4])

        if (c1 == c2 == c3 == c4) & (r1 == r2 == r3 == r4)
            _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                             c1, r1, values, i, n_vars, T(1))
            continue
        end

        # Unwrap corner longitudes to a common branch relative to corner 1 so the
        # bilinear subdivision does not interpolate across the ±180° seam.
        lon1 = Float64(lon_c[i, 1])
        for k in 1:4
            latv[k] = Float64(lat_c[i, k])
            lk = Float64(lon_c[i, k])
            lonv[k] = lk - 360 * round((lk - lon1) / 360)
        end

        compute_subpixels!(points_buf, latv, lonv, n_oversample,
                           lats0_buf, lons0_buf, lats1_buf, lons1_buf)

        for jj in 1:n_oversample, ii in 1:n_oversample
            sub_lat = points_buf[ii, jj, 1]
            sub_lon = points_buf[ii, jj, 2]
            col, row = _spherical_cell(spec, sub_lat, sub_lon)
            _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                             col, row, values, i, n_vars, fac)
        end
    end
    nothing
end

@inline function _sphere_unit_vector(lat, lon)
    sin_lat, cos_lat = sincosd(Float64(lat))
    sin_lon, cos_lon = sincosd(_canonical_longitude(Float64(lon)))
    return (cos_lat * cos_lon, cos_lat * sin_lon, sin_lat)
end

@inline _sphere_dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _sphere_cross(a, b) = (a[2] * b[3] - a[3] * b[2],
                               a[3] * b[1] - a[1] * b[3],
                               a[1] * b[2] - a[2] * b[1])

# Inline the lexicographic comparison of finite, fixed-size sort keys.
@inline function _sphere_isless(a::NTuple{3}, b::NTuple{3})
    return isequal(a[1], b[1]) ?
        (isequal(a[2], b[2]) ? isless(a[3], b[3]) : isless(a[2], b[2])) :
        isless(a[1], b[1])
end

@inline _sphere_isless(a::Tuple{F,NTuple{3,F}}, b::Tuple{F,NTuple{3,F}}) where {F} =
    isequal(a[1], b[1]) ? _sphere_isless(a[2], b[2]) : isless(a[1], b[1])

@inline _sphere_minmax(a, b) = _sphere_isless(b, a) ? (b, a) : (a, b)

@inline function _sort_four(a, b, c, d)
    a, b = _sphere_minmax(a, b)
    c, d = _sphere_minmax(c, d)
    a, c = _sphere_minmax(a, c)
    b, d = _sphere_minmax(b, d)
    b, c = _sphere_minmax(b, c)
    return a, b, c, d
end

@inline function _sphere_ordered_corners(a::NTuple{3,F}, b, c, d) where {F}
    # Summing in a canonical order makes the tangent frame independent of
    # the input permutation, including ties at the angular sort's branch cut.
    a, b, c, d = _sort_four(a, b, c, d)
    centre = a .+ b .+ c .+ d
    norm = _spherical_sqrt(_sphere_dot(centre, centre))
    normal = norm == 0 ? a : centre ./ norm
    axis = abs(normal[3]) > abs(normal[1]) ? (one(F), zero(F), zero(F)) : (zero(F), zero(F), one(F))
    tangent = _sphere_cross(axis, normal)
    tangent_x = tangent ./ _spherical_sqrt(_sphere_dot(tangent, tangent))
    tangent_y = _sphere_cross(normal, tangent_x)
    keyed(point) = (atan(_sphere_dot(point, tangent_y), _sphere_dot(point, tangent_x)), point)
    ordered = _sort_four(keyed(a), keyed(b), keyed(c), keyed(d))
    return map(last, ordered), normal
end

@inline function _sphere_subpoint_edges(corners::NTuple{4,NTuple{3,F}}, j, n) where {F}
    a, b, c, d = corners
    v = (F(j) - F(0.5)) / F(n)
    return (one(F) - v) .* a .+ v .* d, (one(F) - v) .* b .+ v .* c
end

@inline function _sphere_subpoint(left::NTuple{3,F}, right, normal, i, n) where {F}
    u = (F(i) - F(0.5)) / F(n)
    point = (one(F) - u) .* left .+ u .* right
    norm = _spherical_sqrt(_sphere_dot(point, point))
    return iszero(norm) ? normal : point ./ norm
end

@inline function _unit_vector_cell(spec, ::Nothing, x, y, z)
    horizontal = hypot(x, y)
    lat = atand(z, horizontal)
    lon = iszero(horizontal) ? zero(x) : atand(y, x)
    return _spherical_cell(spec, lat, lon)
end

@inline function _unit_vector_cell(spec, locator, x, y, z)
    row, col = divrem(locate(locator, x, y, z) - 1, first(accumulator_shape(spec)))
    return col + 1, row + 1
end

"""
    _accumulate_unit_vector_footprints!(data, m2, weights, compute_std, spec,
                                       lat_c, lon_c, values, n_override, cells)

Sample every footprint on the sphere at the cell centres of its bilinear
parameter grid. The four input corners are in degrees and may be unordered.
Count sub-point hits in each cell and add one weight `k/n²` per cell. This
gives an exact unit weight for a single-cell footprint and the per-footprint
overlap needed for a future sum of squared weights. Reuse the cell-index
buffer `cells`, sized for the largest sample count.
"""
function _accumulate_unit_vector_footprints!(data::AbstractArray{T,3}, m2, weights,
                                             compute_std, spec::AbstractSphericalGridSpec,
                                             lat_c, lon_c, values, n_override, cells) where {T}
    locator = _device_locator(spec)
    cell_spacing = _device_cell_spacing(locator, Float64)
    n_pixels, n_vars = size(values)
    @inbounds for i in 1:n_pixels
        a = _sphere_unit_vector(lat_c[i, 1], lon_c[i, 1])
        b = _sphere_unit_vector(lat_c[i, 2], lon_c[i, 2])
        c = _sphere_unit_vector(lat_c[i, 3], lon_c[i, 3])
        d = _sphere_unit_vector(lat_c[i, 4], lon_c[i, 4])
        corners, normal = _sphere_ordered_corners(a, b, c, d)
        n = n_override === nothing ? _unit_vector_samples(corners, cell_spacing) : n_override
        n_samples = n^2
        for j in 1:n
            left, right = _sphere_subpoint_edges(corners, j, n)
            for k in 1:n
                x, y, z = _sphere_subpoint(left, right, normal, k, n)
                cells[k + (j - 1) * n] = _unit_vector_cell(spec, locator, x, y, z)
            end
        end
        # Sort just this footprint's cell indices in place so repeats are
        # contiguous, with a deterministic order of cell updates.
        sort!(view(cells, 1:n_samples); alg=QuickSort)
        first_hit = 1
        while first_hit <= n_samples
            cell = cells[first_hit]
            next_hit = first_hit + 1
            while next_hit <= n_samples && cells[next_hit] == cell
                next_hit += 1
            end
            col, row = cell
            weight = T((next_hit - first_hit) / n_samples)
            _welford_update!(data, m2, weights, compute_std, col, row, values, i, n_vars, weight)
            first_hit = next_hit
        end
    end
    nothing
end

"""
    accumulate_circular_footprint!(grid_data, grid_std, grid_weights, compute_std,
                                   center_lat_idx, center_lon_idx,
                                   lat_idx, lon_idx, values,
                                   n_pixels, n_vars, n_oversample)

Accumulate observations with circular or near-circular footprints. Symmetric
latitude and longitude radii are inferred from the farthest corner
(`max |corner - center|` along each axis; `≲ 0.1%` edge weight bias for
GOSAT-like viewing geometry), and sampling is performed in `Float64` fractional
grid-index space.

Sub-samples outside `1..n_lon × 1..n_lat` are dropped without renormalizing the
per-cell weight `hits / inside_count`, keeping boundary-cell weights equal to
true fractional footprint overlap so regional boundary means remain
domain-invariant. On global rectangular grids, samples beyond `±180°` longitude
or beyond the poles are dropped rather than wrapped.
"""
function accumulate_circular_footprint!(grid_data::AbstractArray{T,3},
                                        grid_std::Union{Nothing,AbstractArray{T,3}},
                                        grid_weights::AbstractMatrix{T},
                                        compute_std::Bool,
                                        center_lat_idx::AbstractVector,
                                        center_lon_idx::AbstractVector,
                                        lat_idx::AbstractMatrix,
                                        lon_idx::AbstractMatrix,
                                        values::AbstractMatrix,
                                        n_pixels::Int, n_vars::Int,
                                        n_oversample::Int) where {T}
    n_lon, n_lat = size(grid_weights)
    hits = Dict{Tuple{Int,Int},Int}()

    @inbounds for i in 1:n_pixels
        clat, clon, radius_lat, radius_lon =
            _circle_footprint_axes(center_lat_idx[i], center_lon_idx[i],
                                   @view(lat_idx[i, :]), @view(lon_idx[i, :]))

        empty!(hits)
        inside_count = if radius_lat <= eps(Float64) || radius_lon <= eps(Float64)
            0
        else
            _foreach_circle_sample(clat, clon, radius_lat, radius_lon, n_oversample) do lat, lon
                sample_row = floor(Int, lat)
                sample_col = floor(Int, lon)
                if 1 <= sample_col <= n_lon && 1 <= sample_row <= n_lat
                    cell = (sample_col, sample_row)
                    hits[cell] = get(hits, cell, 0) + 1
                end
            end
        end
        if inside_count == 0
            row = floor(Int, clat)
            col = floor(Int, clon)
            if 1 <= col <= n_lon && 1 <= row <= n_lat
                _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                                 col, row, values, i, n_vars, one(T))
            end
        else
            for ((col, row), count) in hits
                weight = T(count / inside_count)
                _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                                 col, row, values, i, n_vars, weight)
            end
        end
    end
    nothing
end

"""
Welford online update for a single grid cell. Adds observation `values[pixel_idx, :]`
with the given `weight` to the running mean at `grid_data[col, row, :]`.
"""
@inline function _welford_update!(grid_data, grid_std, grid_weights, compute_std::Bool,
                                  col::Integer, row::Integer,
                                  values::AbstractMatrix, pixel_idx::Int,
                                  n_vars::Int, weight)
    grid_weights[col, row] += weight
    w = grid_weights[col, row]
    for z in 1:n_vars
        mean_old = grid_data[col, row, z]
        grid_data[col, row, z] = mean_old + weight / w * (values[pixel_idx, z] - mean_old)
        if compute_std
            grid_std[col, row, z] += weight * (values[pixel_idx, z] - mean_old) *
                                     (values[pixel_idx, z] - grid_data[col, row, z])
        end
    end
    nothing
end

"""
    accumulate_batch!(backend, grid_sum, grid_weights,
                      lat_corners, lon_corners, values, n, n_vars)

Accumulate quadrilateral observations into the grid using the KA kernel pipeline.
This function works on KA CPU, CUDA, Metal, and other compatible backends.

The pipeline:
1. Sort corners into CCW order in the input coordinate type (or `Float32` on Metal)
2. Compute `n×n` subpixel grid indices per footprint, mirroring
   [`accumulate_footprint!`](@ref) (`compute_subpixels!` / `floor_indices!`,
   rounding sub-sample positions to the accumulator type `T` before `floor` on
   backends with `Float64` support, and executing in `Float32` on Metal)
3. Scatter-accumulate weighted sums into grid cells (using atomics on GPU)

Sub-pixels outside `1..n_lon × 1..n_lat` are dropped without renormalizing
remaining weights (`1/n²`) so boundary-cell weights reflect true fractional
footprint overlap and regional boundary means remain domain-invariant. On global
rectangular grids, samples beyond `±180°` longitude or beyond the poles are
dropped rather than wrapped.

After processing all files, call [`finalize_mean!`](@ref) to convert sums to means.

This low-level helper accumulates weighted sums in `grid_sum`. Rectangular and spherical
[`grid_l2`](@ref) use a shared reference-shifted accumulator. They replay
the inputs to refine Float32 means and to compute standard deviations when
requested.

# Arguments
- `backend`: KernelAbstractions backend (`CPU()`, `CUDABackend()`, `MetalBackend()`, etc.)
- `grid_sum`: `(n_lon, n_lat, n_vars)` — weighted sum accumulator
- `grid_weights`: `(n_lon, n_lat)` — weight accumulator
- `lat_corners`, `lon_corners`: `(n_fp, 4)` — fractional grid indices of corners
- `values`: `(n_fp, n_vars)` — observation values
- `n`: oversampling factor
- `n_vars`: number of variables being gridded
"""
function accumulate_batch!(backend,
                           grid_sum::AbstractArray{T,3},
                           grid_weights::AbstractMatrix{T},
                           lat_corners::AbstractMatrix,
                           lon_corners::AbstractMatrix,
                           values::AbstractMatrix,
                           n::Int, n_vars::Int) where {T}
    _accumulate_rectangular_batch!(backend, grid_sum, grid_weights, nothing,
        lat_corners, lon_corners, values, n, n_vars, true)
end

"""
    accumulate_circular_batch!(backend, grid_sum, grid_weights,
                               center_lat, center_lon,
                               lat_corners, lon_corners, values, n, n_vars)

Accumulate circular or near-circular footprints using the KA pipeline. This is
the sum-based backend counterpart to [`accumulate_circular_footprint!`](@ref)
and works with KA CPU, CUDA, Metal, and other compatible backends.

Symmetric latitude and longitude radii are inferred from the farthest corner
(`max |corner - center|` along each axis; `≲ 0.1%` edge weight bias for
GOSAT-like viewing geometry). On backends that support double precision
(KA CPU, CUDA), radii and sub-sample positions are promoted to `Float64` inside
the index kernel to mirror `accumulate_circular_footprint!`; on Metal, where
`Float64` is unavailable, they execute in `Float32` and boundary cell
assignments may differ from the sequential path. Each footprint contributes
weight `k / inside_count` to a cell with `k` hits. Sub-samples outside
`1..n_lon × 1..n_lat` are dropped without renormalization, so boundary-cell
weights reflect true fractional footprint overlap and remain domain-invariant.
On global rectangular grids, this also applies beyond `±180°` longitude and
beyond the poles.
"""
function accumulate_circular_batch!(backend,
                                    grid_sum::AbstractArray{T,3},
                                    grid_weights::AbstractMatrix{T},
                                    center_lat::AbstractVector,
                                    center_lon::AbstractVector,
                                    lat_corners::AbstractMatrix,
                                    lon_corners::AbstractMatrix,
                                    values::AbstractMatrix,
                                    n::Int, n_vars::Int) where {T}
    _accumulate_rectangular_batch!(backend, grid_sum, grid_weights, (center_lat, center_lon),
        lat_corners, lon_corners, values, n, n_vars, false)
end

_rect_backend_centers(backend, ::Nothing) = nothing
_rect_backend_centers(backend, centers) = map(a -> _rect_backend_input(backend, a), centers)

function _accumulate_rectangular_batch!(backend, grid_sum, grid_weights, centers,
                                        lat_corners, lon_corners, values, n, n_vars, sort_corners)
    n_fp = size(lat_corners, 1)
    n_fp == 0 && return nothing

    val = _rect_backend_input(backend, values)
    geometry = _rect_batch_geometry(backend, size(grid_weights), eltype(grid_weights),
                                    centers, lat_corners, lon_corners, n, sort_corners)
    (; ix, iy, counts, skip) = geometry
    _scatter_batch!(backend, grid_sum, grid_weights, ix, iy, counts, skip, val, n, n_vars)
    nothing
end

function _rect_batch_geometry(backend, shape, ::Type{P}, centers,
                              lat_corners, lon_corners, n, sort_corners) where {P}
    centers = _rect_backend_centers(backend, centers)
    lat_c = _rect_backend_input(backend, lat_corners)
    lon_c = _rect_backend_input(backend, lon_corners)
    sort_corners && sort_corners_ccw_ka!(backend, lat_c, lon_c)

    # Circles count all inside samples, including those outside the output grid.
    n_fp = size(lat_corners, 1)
    ix = KernelAbstractions.allocate(backend, Int32, (n_fp, n * n))
    iy = KernelAbstractions.allocate(backend, Int32, (n_fp, n * n))
    skip = KernelAbstractions.allocate(backend, Int32, (n_fp,))
    counts = _rect_batch_indices!(backend, shape, P, ix, iy, skip, centers, lat_c, lon_c, n)
    (; ix, iy, counts, skip)
end

function _rect_batch_indices!(backend, shape, ::Type{P}, ix, iy, skip, ::Nothing, lat, lon, n) where {P}
    n_lon, n_lat = shape
    _compute_footprint_indices_ka!(backend, ix, iy, skip, lat, lon, n,
        n_lon, n_lat, _ka_supported_float(backend, P))
    nothing
end

function _rect_batch_indices!(backend, shape, ::Type{P}, ix, iy, skip, centers, lat, lon, n) where {P}
    n_lon, n_lat = shape
    counts = KernelAbstractions.allocate(backend, Int32, (size(lat, 1),))
    compute_circular_footprint_indices_ka!(backend, ix, iy, counts, skip,
        centers..., lat, lon, n; n_lon, n_lat)
    counts
end

"""
    finalize_mean!(grid_data, grid_weights)

Convert weighted-sum accumulator to mean: `mean = sum / weight`.
Cells with near-zero weight are left as zero.
"""
function finalize_mean!(grid_data::AbstractArray{T,3}, grid_weights::AbstractMatrix{T}) where {T}
    n_cells = length(grid_weights)
    @inbounds for cell in 1:n_cells
        _finalize_mean_cell!(grid_data, grid_weights, cell, n_cells, size(grid_data, 3), eps(T))
    end
    nothing
end

@inline function _finalize_mean_cell!(sums, weights, cell, stride, n_vars, min_weight)
    @inbounds w = weights[cell]
    if w > min_weight
        @inbounds for v in 1:n_vars
            sums[cell + (v - 1) * stride] /= w
        end
    end
    return nothing
end

"""
    accumulate_center!(grid_data, grid_weights, lat, lon, values, grid_spec)

Accumulate observations using center coordinates only (no footprint bounds).
Used for MODIS-style gridding where each pixel maps to exactly one grid cell.
Simple sum accumulation — divide by `grid_weights` afterward to get the mean.
"""
function accumulate_center!(grid_data::AbstractArray{T,3},
                            grid_weights::AbstractMatrix{T},
                            lat_idx::AbstractVector{<:Integer},
                            lon_idx::AbstractVector{<:Integer},
                            values::AbstractMatrix{T},
                            n_pixels::Int, n_vars::Int) where {T}
    _accumulate_center!(grid_data, grid_weights, lat_idx, lon_idx, values, nothing, n_pixels, n_vars)
end

# Mask-aware variant: skip invalid rows in place rather than allocating
# compacted copies of `lat_idx`/`lon_idx`/`values` upstream.
function accumulate_center!(grid_data::AbstractArray{T,3},
                            grid_weights::AbstractMatrix{T},
                            lat_idx::AbstractVector{<:Integer},
                            lon_idx::AbstractVector{<:Integer},
                            values::AbstractMatrix{T},
                            valid::BitVector,
                            n_pixels::Int, n_vars::Int) where {T}
    _accumulate_center!(grid_data, grid_weights, lat_idx, lon_idx, values, valid, n_pixels, n_vars)
end

@inline _valid_center_row(::Nothing, i) = true
@inline _valid_center_row(valid, i) = @inbounds valid[i]

function _accumulate_center!(grid_data, grid_weights, lat_idx, lon_idx, values, valid, n_pixels, n_vars)
    T = eltype(grid_weights)
    @inbounds for i in 1:n_pixels
        _valid_center_row(valid, i) || continue
        col = lon_idx[i]
        row = lat_idx[i]
        grid_weights[col, row] += one(T)
        for z in 1:n_vars
            grid_data[col, row, z] += values[i, z]
        end
    end
    nothing
end

"""
    finalize_std!(grid_std, grid_weights)

Convert Welford M2 accumulator to standard deviation: `std = sqrt(M2 / weight)`.
"""
function finalize_std!(grid_std::AbstractArray{T,3}, grid_weights::AbstractMatrix{T}) where {T}
    nlon, nlat, nvars = size(grid_std)
    @inbounds for j in 1:nlat, i in 1:nlon
        w = grid_weights[i, j]
        if w > eps(T)
            for z in 1:nvars
                grid_std[i, j, z] = sqrt(grid_std[i, j, z] / w)
            end
        end
    end
    nothing
end
