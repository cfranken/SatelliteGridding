"""
    sort_corners_ccw(lat_corners, lon_corners) -> (sorted_lat, sorted_lon)

Sort 4 footprint corners into counter-clockwise order by computing the angle from
the centroid to each corner. This ensures edges 1→2 and 4→3 are always opposite
edges, which is required by `compute_subpixels!`.

Different satellite products store corners in different orders (e.g., TROPOMI uses
CCW order while OCO-2 uses a different convention). This function normalizes them.
"""
function sort_corners_ccw(lat_corners::AbstractVector, lon_corners::AbstractVector)
    clat = (lat_corners[1] + lat_corners[2] + lat_corners[3] + lat_corners[4]) / 4
    clon = (lon_corners[1] + lon_corners[2] + lon_corners[3] + lon_corners[4]) / 4
    angles = (atan(lat_corners[1] - clat, lon_corners[1] - clon),
              atan(lat_corners[2] - clat, lon_corners[2] - clon),
              atan(lat_corners[3] - clat, lon_corners[3] - clon),
              atan(lat_corners[4] - clat, lon_corners[4] - clon))
    order = sortperm(collect(angles))
    return lat_corners[order], lon_corners[order]
end

"""
    sort_corners_ccw!(lat_mat, lon_mat)

Sort corners in-place for all soundings in N×4 matrices.
"""
function sort_corners_ccw!(lat_mat::AbstractMatrix, lon_mat::AbstractMatrix)
    n = size(lat_mat, 1)
    @inbounds for i in 1:n
        clat = (lat_mat[i,1] + lat_mat[i,2] + lat_mat[i,3] + lat_mat[i,4]) / 4
        clon = (lon_mat[i,1] + lon_mat[i,2] + lon_mat[i,3] + lon_mat[i,4]) / 4
        a1 = atan(lat_mat[i,1] - clat, lon_mat[i,1] - clon)
        a2 = atan(lat_mat[i,2] - clat, lon_mat[i,2] - clon)
        a3 = atan(lat_mat[i,3] - clat, lon_mat[i,3] - clon)
        a4 = atan(lat_mat[i,4] - clat, lon_mat[i,4] - clon)

        # Simple sorting network for 4 elements
        la = (lat_mat[i,1], lat_mat[i,2], lat_mat[i,3], lat_mat[i,4])
        lo = (lon_mat[i,1], lon_mat[i,2], lon_mat[i,3], lon_mat[i,4])
        order = sortperm([a1, a2, a3, a4])

        lat_mat[i,1] = la[order[1]]; lon_mat[i,1] = lo[order[1]]
        lat_mat[i,2] = la[order[2]]; lon_mat[i,2] = lo[order[2]]
        lat_mat[i,3] = la[order[3]]; lon_mat[i,3] = lo[order[3]]
        lat_mat[i,4] = la[order[4]]; lon_mat[i,4] = lo[order[4]]
    end
    nothing
end

"""
    compute_n_oversample(footprint_extent, cell_size; n_min=2, n_max=20)

Estimate an appropriate oversampling factor `n` based on the ratio of footprint size
to grid cell size. Aims for ~3 sub-pixels per grid cell width.

Returns an integer in `[n_min, n_max]`.
"""
function compute_n_oversample(footprint_extent::Real, cell_size::Real;
                              n_min::Int=2, n_max::Int=20)
    ratio = footprint_extent / cell_size
    clamp(ceil(Int, ratio * 3), n_min, n_max)
end

"""
    subdivide_line!(points, j, lat1, lon1, lat2, lon2, n)

Subdivide the line segment from `(lat1, lon1)` to `(lat2, lon2)` into `n`
equally-spaced interior points. Results are stored in `points[j, 1:n, 1:2]`.

The points are placed at the centers of `n` equal sub-segments along the line.
"""
function subdivide_line!(points::AbstractArray{T,3}, j::Int,
                         lat1::Real, lon1::Real, lat2::Real, lon2::Real,
                         n::Int) where {T}
    dlat = (lat2 - lat1) / (2n)
    dlon = (lon2 - lon1) / (2n)
    start_lat = lat1 + dlat
    start_lon = lon1 + dlon
    @inbounds for i in 1:n
        points[j, i, 1] = start_lat + 2 * (i - 1) * dlat
        points[j, i, 2] = start_lon + 2 * (i - 1) * dlon
    end
    nothing
end

"""
    subdivide_baseline!(lats, lons, lat1, lon1, lat2, lon2, n)

Subdivide a baseline (edge) into `n` interior points, storing lat and lon separately.
Used as the first step in `compute_subpixels!`.
"""
function subdivide_baseline!(lats::AbstractVector, lons::AbstractVector,
                             lat1::Real, lon1::Real, lat2::Real, lon2::Real,
                             n::Int)
    dlat = (lat2 - lat1) / (2n)
    dlon = (lon2 - lon1) / (2n)
    start_lat = lat1 + dlat
    start_lon = lon1 + dlon
    @inbounds for i in 1:n
        lats[i] = start_lat + 2 * (i - 1) * dlat
        lons[i] = start_lon + 2 * (i - 1) * dlon
    end
    nothing
end

"""
    compute_subpixels!(points, vert_lat, vert_lon, n, lats_0, lons_0, lats_1, lons_1)

Fill the `n×n×2` array `points` with sub-pixel positions inside the quadrilateral
defined by 4 corner coordinates (`vert_lat[1:4]`, `vert_lon[1:4]`).

The algorithm:
1. Subdivide edge 1→2 into `n` baseline points
2. Subdivide edge 4→3 into `n` baseline points
3. For each pair of corresponding baseline points, subdivide the connecting
   line into `n` interior points

This produces `n²` sub-pixels that fill the footprint quadrilateral.
Temporary buffers `lats_0`, `lons_0`, `lats_1`, `lons_1` (each length `n`) must be provided.
"""
function compute_subpixels!(points::AbstractArray{T,3},
                            vert_lat::AbstractVector, vert_lon::AbstractVector,
                            n::Int,
                            lats_0::AbstractVector, lons_0::AbstractVector,
                            lats_1::AbstractVector, lons_1::AbstractVector) where {T}
    # Subdivide the two opposite edges
    subdivide_baseline!(lats_0, lons_0, vert_lat[1], vert_lon[1], vert_lat[2], vert_lon[2], n)
    subdivide_baseline!(lats_1, lons_1, vert_lat[4], vert_lon[4], vert_lat[3], vert_lon[3], n)

    # Connect corresponding points across the two edges
    for i in 1:n
        subdivide_line!(points, i, lats_0[i], lons_0[i], lats_1[i], lons_1[i], n)
    end
    nothing
end

"""
    floor_indices!(ix, iy, points)

Convert the floating-point sub-pixel positions in `points` to integer grid cell indices
via `floor`. Results stored in `ix` and `iy` (each length `n²`).
"""
function floor_indices!(ix::AbstractVector{<:Integer}, iy::AbstractVector{<:Integer},
                        points::AbstractArray)
    points_x = @view points[:, :, 1]
    points_y = @view points[:, :, 2]
    @inbounds for ic in LinearIndices(points_x)
        ix[ic] = floor(Int32, points_x[ic])
        iy[ic] = floor(Int32, points_y[ic])
    end
    nothing
end

# --- KernelAbstractions kernels for batch processing (CPU + GPU) ---

"""
    sort_corners_ccw_ka!(backend, lat_mat, lon_mat)

KA kernel that sorts footprint corners into CCW order for all soundings in parallel.
Uses a 5-comparator sorting network (optimal for 4 elements, no allocations).

- `lat_mat`, `lon_mat`: N×4 matrices of corner coordinates
"""
function sort_corners_ccw_ka!(backend, lat_mat::AbstractMatrix, lon_mat::AbstractMatrix)
    n = size(lat_mat, 1)
    kernel! = _sort_ccw_kernel!(backend)
    kernel!(lat_mat, lon_mat, ndrange=n)
    KernelAbstractions.synchronize(backend)
    nothing
end

@kernel function _sort_ccw_kernel!(lat_mat, lon_mat)
    i = @index(Global)

    @inbounds begin
        # Compute centroid
        clat = (lat_mat[i,1] + lat_mat[i,2] + lat_mat[i,3] + lat_mat[i,4]) / 4
        clon = (lon_mat[i,1] + lon_mat[i,2] + lon_mat[i,3] + lon_mat[i,4]) / 4

        # Compute angles
        a1 = atan(lat_mat[i,1] - clat, lon_mat[i,1] - clon)
        a2 = atan(lat_mat[i,2] - clat, lon_mat[i,2] - clon)
        a3 = atan(lat_mat[i,3] - clat, lon_mat[i,3] - clon)
        a4 = atan(lat_mat[i,4] - clat, lon_mat[i,4] - clon)

        # Load corners
        la1, la2, la3, la4 = lat_mat[i,1], lat_mat[i,2], lat_mat[i,3], lat_mat[i,4]
        lo1, lo2, lo3, lo4 = lon_mat[i,1], lon_mat[i,2], lon_mat[i,3], lon_mat[i,4]

        # 5-comparator sorting network for 4 elements
        if a1 > a2
            a1, a2 = a2, a1; la1, la2 = la2, la1; lo1, lo2 = lo2, lo1
        end
        if a3 > a4
            a3, a4 = a4, a3; la3, la4 = la4, la3; lo3, lo4 = lo4, lo3
        end
        if a1 > a3
            a1, a3 = a3, a1; la1, la3 = la3, la1; lo1, lo3 = lo3, lo1
        end
        if a2 > a4
            a2, a4 = a4, a2; la2, la4 = la4, la2; lo2, lo4 = lo4, lo2
        end
        if a2 > a3
            a2, a3 = a3, a2; la2, la3 = la3, la2; lo2, lo3 = lo3, lo2
        end

        # Write back sorted corners
        lat_mat[i,1] = la1; lat_mat[i,2] = la2; lat_mat[i,3] = la3; lat_mat[i,4] = la4
        lon_mat[i,1] = lo1; lon_mat[i,2] = lo2; lon_mat[i,3] = lo3; lon_mat[i,4] = lo4
    end
end

"""
    compute_footprint_indices_ka!(backend, ix_out, iy_out, skip_flag,
                                  lat_corners, lon_corners, n)

KA kernel that computes n×n sub-pixel grid indices for each footprint.
Also sets `skip_flag[fp] = true` if the footprint's longitude extent >= n
(too wide to oversample meaningfully).

For footprints where all 4 corners fall in the same cell, the kernel
stores the single cell index in position 1 and marks all other positions
with ix=-1 (sentinel). The accumulation kernel handles these fast-path
footprints by checking `skip_flag`.

- `ix_out`: (n_fp, n²) Int32 — lat grid cell indices
- `iy_out`: (n_fp, n²) Int32 — lon grid cell indices
- `skip_flag`: (n_fp,) Int32 — 0=oversample, 1=fast-path (single cell), 2=skip (too wide)
- `lat_corners`, `lon_corners`: (n_fp, 4) — fractional grid indices of corners
- `n`: oversampling factor
"""
function compute_footprint_indices_ka!(backend,
                                       ix_out::AbstractMatrix{Int32},
                                       iy_out::AbstractMatrix{Int32},
                                       skip_flag::AbstractVector{Int32},
                                       lat_corners::AbstractMatrix,
                                       lon_corners::AbstractMatrix,
                                       n::Int)
    n_fp = size(lat_corners, 1)
    kernel! = _compute_indices_kernel!(backend)
    kernel!(ix_out, iy_out, skip_flag, lat_corners, lon_corners, Int32(n), ndrange=n_fp)
    KernelAbstractions.synchronize(backend)
    nothing
end

@kernel function _compute_indices_kernel!(ix_out, iy_out, skip_flag,
                                          lat_corners, lon_corners, n::Int32)
    fp = @index(Global)

    @inbounds begin
        # Compute grid-cell bounding box for this footprint
        il1 = floor(Int32, lat_corners[fp, 1])
        il2 = floor(Int32, lat_corners[fp, 2])
        il3 = floor(Int32, lat_corners[fp, 3])
        il4 = floor(Int32, lat_corners[fp, 4])
        io1 = floor(Int32, lon_corners[fp, 1])
        io2 = floor(Int32, lon_corners[fp, 2])
        io3 = floor(Int32, lon_corners[fp, 3])
        io4 = floor(Int32, lon_corners[fp, 4])

        min_lat = min(il1, il2, il3, il4)
        max_lat = max(il1, il2, il3, il4)
        min_lon = min(io1, io2, io3, io4)
        max_lon = max(io1, io2, io3, io4)

        dim_lat = max_lat - min_lat
        dim_lon = max_lon - min_lon

        nsq = n * n

        if (dim_lat == Int32(0)) & (dim_lon == Int32(0))
            # Fast path: all corners in one cell
            skip_flag[fp] = Int32(1)
            ix_out[fp, 1] = io1
            iy_out[fp, 1] = il1
            for k in Int32(2):nsq
                ix_out[fp, k] = Int32(-1)
                iy_out[fp, k] = Int32(-1)
            end
        elseif dim_lon >= n
            # Too wide: skip this footprint
            skip_flag[fp] = Int32(2)
            for k in Int32(1):nsq
                ix_out[fp, k] = Int32(-1)
                iy_out[fp, k] = Int32(-1)
            end
        else
            # Oversample: compute n×n subpixel indices
            skip_flag[fp] = Int32(0)

            # Subdivide edge 1→2 and edge 4→3, then connect across
            inv2n = one(eltype(lat_corners)) / (2 * n)
            for i in Int32(1):n
                # Points along edge 1→2
                t0 = (2 * i - Int32(1)) * inv2n
                lat0 = lat_corners[fp, 1] + t0 * (lat_corners[fp, 2] - lat_corners[fp, 1])
                lon0 = lon_corners[fp, 1] + t0 * (lon_corners[fp, 2] - lon_corners[fp, 1])
                # Points along edge 4→3
                lat1 = lat_corners[fp, 4] + t0 * (lat_corners[fp, 3] - lat_corners[fp, 4])
                lon1 = lon_corners[fp, 4] + t0 * (lon_corners[fp, 3] - lon_corners[fp, 4])

                for j in Int32(1):n
                    t1 = (2 * j - Int32(1)) * inv2n
                    sub_lat = lat0 + t1 * (lat1 - lat0)
                    sub_lon = lon0 + t1 * (lon1 - lon0)
                    k = (i - Int32(1)) * n + j
                    ix_out[fp, k] = floor(Int32, sub_lon)
                    iy_out[fp, k] = floor(Int32, sub_lat)
                end
            end
        end
    end
end

"""
    scatter_accumulate!(grid_sum, grid_weights, ix, iy, skip_flag, values, n, n_vars)

Sequential scatter-accumulate of weighted values into grid cells.
Uses precomputed grid indices from `compute_footprint_indices_ka!`.

For each footprint:
- `skip_flag == 1` (fast path): add weight=1, value to single cell
- `skip_flag == 0` (oversample): add weight=1/n² to each subpixel's cell
- `skip_flag == 2` (skip): do nothing

This function is sequential (no atomics needed) and works with any array type.
For GPU backends, use the `@kernel` version with `@atomic` instead.

- `grid_sum`: (n_lon, n_lat, n_vars) — weighted sum accumulator
- `grid_weights`: (n_lon, n_lat) — weight accumulator
- `ix`, `iy`: (n_fp, n²) — precomputed grid indices
- `skip_flag`: (n_fp,) — 0=oversample, 1=fast-path, 2=skip
- `values`: (n_fp, n_vars) — observation values
"""
function scatter_accumulate!(grid_sum::AbstractArray{T,3},
                             grid_weights::AbstractMatrix{T},
                             ix::AbstractMatrix{<:Integer},
                             iy::AbstractMatrix{<:Integer},
                             skip_flag::AbstractVector{<:Integer},
                             values::AbstractMatrix,
                             n::Int, n_vars::Int) where {T}
    n_fp = size(ix, 1)
    fac = T(1 / n^2)
    nsq = n * n

    @inbounds for fp in 1:n_fp
        flag = skip_flag[fp]

        if flag == 1
            # Fast path: single cell, weight = 1
            col = Int(ix[fp, 1])
            row = Int(iy[fp, 1])
            grid_weights[col, row] += one(T)
            for z in 1:n_vars
                grid_sum[col, row, z] += values[fp, z]
            end

        elseif flag == 0
            # Oversample path: n² subpixels, weight = 1/n²
            for k in 1:nsq
                col = Int(ix[fp, k])
                row = Int(iy[fp, k])
                grid_weights[col, row] += fac
                for z in 1:n_vars
                    grid_sum[col, row, z] += fac * values[fp, z]
                end
            end
        end
        # flag == 2: skip (footprint too wide)
    end
    nothing
end

"""
    scatter_accumulate_ka!(backend, grid_sum, grid_weights,
                           ix, iy, skip_flag, values, n, n_vars)

KA kernel version of scatter-accumulate using `@atomic` for thread safety.
Required for GPU backends where multiple threads write to the same grid cells.
For CPU backends, use `scatter_accumulate!` (sequential) instead.
"""
function scatter_accumulate_ka!(backend,
                                grid_sum::AbstractArray{T,3},
                                grid_weights::AbstractMatrix{T},
                                ix::AbstractMatrix{Int32},
                                iy::AbstractMatrix{Int32},
                                skip_flag::AbstractVector{Int32},
                                values::AbstractMatrix,
                                n::Int, n_vars::Int) where {T}
    n_fp = size(ix, 1)
    kernel! = _scatter_kernel!(backend)
    kernel!(grid_sum, grid_weights, ix, iy, skip_flag, values,
            Int32(n), Int32(n_vars), ndrange=n_fp)
    KernelAbstractions.synchronize(backend)
    nothing
end

@kernel function _scatter_kernel!(grid_sum, grid_weights, ix, iy, skip_flag,
                                  values, n::Int32, n_vars::Int32)
    fp = @index(Global)

    @inbounds begin
        flag = skip_flag[fp]

        if flag == Int32(1)
            # Fast path: single cell, weight = 1
            col = ix[fp, 1]
            row = iy[fp, 1]
            @atomic grid_weights[col, row] += one(eltype(grid_weights))
            for z in Int32(1):n_vars
                @atomic grid_sum[col, row, z] += values[fp, z]
            end

        elseif flag == Int32(0)
            # Oversample path: n² subpixels, weight = 1/n²
            fac = one(eltype(grid_weights)) / (n * n)
            nsq = n * n
            for k in Int32(1):nsq
                col = ix[fp, k]
                row = iy[fp, k]
                @atomic grid_weights[col, row] += fac
                for z in Int32(1):n_vars
                    @atomic grid_sum[col, row, z] += fac * values[fp, z]
                end
            end
        end
        # flag == 2: skip (footprint too wide)
    end
end
