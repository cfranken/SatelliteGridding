"""
    accumulate_footprint!(grid_data, grid_std, grid_weights, compute_std,
                          lat_idx, lon_idx, values, n_pixels, n_vars, n_oversample,
                          grid_spec, points_buf, ix_buf, iy_buf,
                          lats0_buf, lons0_buf, lats1_buf, lons1_buf)

Accumulate a batch of satellite observations into the output grid using footprint-aware
oversampling and Welford's online averaging algorithm.

For each pixel:
- **Fast path**: If all 4 corner coordinates fall into the same grid cell, the full
  observation value is assigned to that cell (weight=1).
- **Oversample path**: If the footprint spans multiple cells, it is subdivided into
  `n_oversample × n_oversample` sub-pixels. Each sub-pixel contributes weight `1/n²`
  to whichever grid cell it falls in.

The Welford algorithm maintains a running mean and sum-of-squared-deviations (M2),
enabling numerically stable computation of mean and variance in a single pass.

# Arguments
- `grid_data`: Output array `(n_lon, n_lat, n_vars)` — running mean
- `grid_std`: Output array `(n_lon, n_lat, n_vars)` — running M2 (sum of squared deviations)
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
                               grid_std::AbstractArray{T,3},
                               grid_weights::AbstractMatrix{T},
                               compute_std::Bool,
                               lat_idx::AbstractMatrix, lon_idx::AbstractMatrix,
                               values::AbstractMatrix,
                               n_pixels::Int, n_vars::Int, n_oversample::Int,
                               points_buf::AbstractArray{T,3},
                               ix_buf::AbstractVector{<:Integer},
                               iy_buf::AbstractVector{<:Integer},
                               lats0_buf::AbstractVector,
                               lons0_buf::AbstractVector,
                               lats1_buf::AbstractVector,
                               lons1_buf::AbstractVector) where {T}
    # Compute grid-cell indices for all 4 corners
    ilon = floor.(Int32, lon_idx)
    ilat = floor.(Int32, lat_idx)
    min_lat = minimum(Int32, floor.(lat_idx), dims=2)
    max_lat = maximum(Int32, floor.(lat_idx), dims=2)
    min_lon = minimum(Int32, floor.(lon_idx), dims=2)
    max_lon = maximum(Int32, floor.(lon_idx), dims=2)
    dim_lat = max_lat .- min_lat
    dim_lon = max_lon .- min_lon
    dist_lon = max_lon .- min_lon

    fac = T(1 / n_oversample^2)

    @inbounds for i in 1:n_pixels
        if (dim_lat[i] == 0) & (dim_lon[i] == 0)
            # Fast path: all corners in one cell
            _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                             ilon[i, 1], ilat[i, 1], values, i, n_vars, T(1))
        elseif dist_lon[i] < n_oversample
            # Oversample path: footprint spans multiple cells
            compute_subpixels!(points_buf, lat_idx[i, :], lon_idx[i, :], n_oversample,
                               lats0_buf, lons0_buf, lats1_buf, lons1_buf)
            floor_indices!(ix_buf, iy_buf, points_buf)

            @inbounds for j in eachindex(ix_buf)
                _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                                 iy_buf[j], ix_buf[j], values, i, n_vars, fac)
            end
        end
    end
    nothing
end

"""
    accumulate_circular_footprint!(grid_data, grid_std, grid_weights, compute_std,
                                   center_lat_idx, center_lon_idx,
                                   lat_idx, lon_idx, values,
                                   n_pixels, n_vars, n_oversample)

Accumulate observations with circular or near-circular footprints. The center
coordinate gives the footprint center and the four coordinates define its
bounding box or edge points. Sampling is done in fractional grid-index space.
"""
function accumulate_circular_footprint!(grid_data::AbstractArray{T,3},
                                        grid_std::AbstractArray{T,3},
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

    @inbounds for i in 1:n_pixels
        clat, clon, radius_lat, radius_lon =
            _circle_footprint_axes(center_lat_idx[i], center_lon_idx[i],
                                   @view(lat_idx[i, :]), @view(lon_idx[i, :]))

        if radius_lat <= eps(Float64) || radius_lon <= eps(Float64)
            row = floor(Int, clat)
            col = floor(Int, clon)
            if 1 <= col <= n_lon && 1 <= row <= n_lat
                _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                                 col, row, values, i, n_vars, T(1))
            end
            continue
        end

        inside_count = 0
        for iy in 1:n_oversample, ix in 1:n_oversample
            lat = clat - radius_lat + (iy - 0.5) * (2radius_lat / n_oversample)
            lon = clon - radius_lon + (ix - 0.5) * (2radius_lon / n_oversample)
            if _inside_bounded_circle(lat, lon, clat, clon, radius_lat, radius_lon)
                inside_count += 1
            end
        end

        if inside_count == 0
            row = floor(Int, clat)
            col = floor(Int, clon)
            if 1 <= col <= n_lon && 1 <= row <= n_lat
                _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                                 col, row, values, i, n_vars, T(1))
            end
            continue
        end

        weight = T(1 / inside_count)
        for iy in 1:n_oversample, ix in 1:n_oversample
            lat = clat - radius_lat + (iy - 0.5) * (2radius_lat / n_oversample)
            lon = clon - radius_lon + (ix - 0.5) * (2radius_lon / n_oversample)
            if _inside_bounded_circle(lat, lon, clat, clon, radius_lat, radius_lon)
                row = floor(Int, lat)
                col = floor(Int, lon)
                if 1 <= col <= n_lon && 1 <= row <= n_lat
                    _welford_update!(grid_data, grid_std, grid_weights, compute_std,
                                     col, row, values, i, n_vars, weight)
                end
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
                      lat_corners, lon_corners, values, n;
                      compute_std=false, grid_std=nothing, grid_mean=nothing)

Accumulate observations into the grid using the KA kernel pipeline.
This function works identically on CPU and GPU backends.

The pipeline:
1. Sort corners into CCW order
2. Compute n×n subpixel grid indices per footprint
3. Scatter-accumulate weighted sums into grid cells (using atomics)

After processing all files, call `finalize_mean!` to convert sums to means.

For `compute_std=true`, a two-pass approach is used:
1. First call `accumulate_batch!` without std to accumulate sums
2. Call `finalize_mean!` to get the mean
3. Call `accumulate_std_pass!` with the computed mean to accumulate variance

# Arguments
- `backend`: KernelAbstractions backend (`CPU()`, `CUDABackend()`, `MetalBackend()`, etc.)
- `grid_sum`: (n_lon, n_lat, n_vars) — weighted sum accumulator
- `grid_weights`: (n_lon, n_lat) — weight accumulator
- `lat_corners`, `lon_corners`: (n_fp, 4) — fractional grid indices of corners
- `values`: (n_fp, n_vars) — observation values
- `n`: oversampling factor
"""
function accumulate_batch!(backend,
                           grid_sum::AbstractArray{T,3},
                           grid_weights::AbstractMatrix{T},
                           lat_corners::AbstractMatrix,
                           lon_corners::AbstractMatrix,
                           values::AbstractMatrix,
                           n::Int, n_vars::Int) where {T}
    n_fp = size(lat_corners, 1)
    n_fp == 0 && return nothing

    # Transfer inputs to backend (no-op for CPU — arrays are already there)
    if backend isa CPU
        lat_c = lat_corners
        lon_c = lon_corners
        val = values
    else
        lat_c = KernelAbstractions.allocate(backend, T, size(lat_corners))
        copyto!(lat_c, lat_corners)
        lon_c = KernelAbstractions.allocate(backend, T, size(lon_corners))
        copyto!(lon_c, lon_corners)
        val = KernelAbstractions.allocate(backend, eltype(values), size(values))
        copyto!(val, values)
    end

    # Step 1: Sort corners into CCW order (KA kernel — parallel)
    sort_corners_ccw_ka!(backend, lat_c, lon_c)

    # Step 2: Compute subpixel indices (KA kernel — parallel)
    nsq = n * n
    ix = KernelAbstractions.zeros(backend, Int32, n_fp, nsq)
    iy = KernelAbstractions.zeros(backend, Int32, n_fp, nsq)
    skip = KernelAbstractions.zeros(backend, Int32, n_fp)
    compute_footprint_indices_ka!(backend, ix, iy, skip, lat_c, lon_c, n)

    # Step 3: Scatter-accumulate
    # CPU: sequential scatter (no atomics needed, better cache locality)
    # GPU: use @atomic kernel for parallel scatter
    if backend isa CPU
        scatter_accumulate!(grid_sum, grid_weights, ix, iy, skip, val, n, n_vars)
    else
        scatter_accumulate_ka!(backend, grid_sum, grid_weights, ix, iy, skip, val, n, n_vars)
    end

    nothing
end

"""
    accumulate_circular_batch!(backend, grid_sum, grid_weights,
                               center_lat, center_lon,
                               lat_corners, lon_corners, values, n, n_vars)

Accumulate circular or near-circular footprints using the KA pipeline. This is
the sum-based backend counterpart to `accumulate_circular_footprint!` and works
with KA CPU, CUDA, Metal, and other compatible backends.
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
    n_fp = size(lat_corners, 1)
    n_fp == 0 && return nothing

    if backend isa CPU
        clat = center_lat
        clon = center_lon
        lat_c = lat_corners
        lon_c = lon_corners
        val = values
    else
        clat = KernelAbstractions.allocate(backend, T, size(center_lat))
        copyto!(clat, center_lat)
        clon = KernelAbstractions.allocate(backend, T, size(center_lon))
        copyto!(clon, center_lon)
        lat_c = KernelAbstractions.allocate(backend, T, size(lat_corners))
        copyto!(lat_c, lat_corners)
        lon_c = KernelAbstractions.allocate(backend, T, size(lon_corners))
        copyto!(lon_c, lon_corners)
        val = KernelAbstractions.allocate(backend, eltype(values), size(values))
        copyto!(val, values)
    end

    nsq = n * n
    ix = KernelAbstractions.zeros(backend, Int32, n_fp, nsq)
    iy = KernelAbstractions.zeros(backend, Int32, n_fp, nsq)
    inside_count = KernelAbstractions.zeros(backend, Int32, n_fp)
    skip = KernelAbstractions.zeros(backend, Int32, n_fp)
    compute_circular_footprint_indices_ka!(backend, ix, iy, inside_count, skip,
                                           clat, clon, lat_c, lon_c, n)

    if backend isa CPU
        scatter_accumulate_circular!(grid_sum, grid_weights, ix, iy, inside_count,
                                     skip, val, n, n_vars)
    else
        scatter_accumulate_circular_ka!(backend, grid_sum, grid_weights, ix, iy,
                                        inside_count, skip, val, n, n_vars)
    end

    nothing
end

"""
    finalize_mean!(grid_data, grid_weights)

Convert weighted-sum accumulator to mean: `mean = sum / weight`.
Cells with near-zero weight are left as zero.
"""
function finalize_mean!(grid_data::AbstractArray{T,3}, grid_weights::AbstractMatrix{T}) where {T}
    nlon, nlat, nvars = size(grid_data)
    @inbounds for j in 1:nlat, i in 1:nlon
        w = grid_weights[i, j]
        if w > eps(T)
            for z in 1:nvars
                grid_data[i, j, z] /= w
            end
        end
    end
    nothing
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
    @inbounds for i in 1:n_pixels
        grid_weights[lon_idx[i], lat_idx[i]] += one(T)
        for z in 1:n_vars
            grid_data[lon_idx[i], lat_idx[i], z] += values[i, z]
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
