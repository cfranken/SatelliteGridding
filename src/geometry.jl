"""
    polygon_area(points)

Compute planar polygon area with the shoelace formula. Points are `(x, y)`
tuples in grid-index coordinates.
"""
function polygon_area(points::AbstractVector)
    n = length(points)
    n < 3 && return 0.0

    area2 = 0.0
    @inbounds for i in 1:n
        x1, y1 = points[i]
        x2, y2 = points[i == n ? 1 : i + 1]
        area2 += x1 * y2 - x2 * y1
    end
    abs(area2) / 2
end

function _line_intersection_x(p1, p2, x)
    x1, y1 = p1
    x2, y2 = p2
    denom = x2 - x1
    abs(denom) <= eps(Float64) && return (x, y2)
    t = (x - x1) / denom
    (x, y1 + t * (y2 - y1))
end

function _line_intersection_y(p1, p2, y)
    x1, y1 = p1
    x2, y2 = p2
    denom = y2 - y1
    abs(denom) <= eps(Float64) && return (x2, y)
    t = (y - y1) / denom
    (x1 + t * (x2 - x1), y)
end

function _clip_polygon(poly, inside, intersection)
    isempty(poly) && return poly

    out = Tuple{Float64,Float64}[]
    prev = poly[end]
    prev_inside = inside(prev)

    for curr in poly
        curr_inside = inside(curr)
        if curr_inside
            !prev_inside && push!(out, intersection(prev, curr))
            push!(out, curr)
        elseif prev_inside
            push!(out, intersection(prev, curr))
        end
        prev = curr
        prev_inside = curr_inside
    end

    out
end

function _clip_polygon_to_cell(poly, col::Integer, row::Integer)
    xmin = Float64(col)
    xmax = Float64(col + 1)
    ymin = Float64(row)
    ymax = Float64(row + 1)

    clipped = _clip_polygon(poly, p -> p[1] >= xmin,
                            (p1, p2) -> _line_intersection_x(p1, p2, xmin))
    clipped = _clip_polygon(clipped, p -> p[1] <= xmax,
                            (p1, p2) -> _line_intersection_x(p1, p2, xmax))
    clipped = _clip_polygon(clipped, p -> p[2] >= ymin,
                            (p1, p2) -> _line_intersection_y(p1, p2, ymin))
    _clip_polygon(clipped, p -> p[2] <= ymax,
                  (p1, p2) -> _line_intersection_y(p1, p2, ymax))
end

function _footprint_polygon(lat_corners, lon_corners)
    lat_sorted, lon_sorted = sort_corners_ccw(collect(lat_corners), collect(lon_corners))
    [(Float64(lon_sorted[i]), Float64(lat_sorted[i])) for i in eachindex(lat_sorted)]
end

"""
    exact_footprint_weights(lat_corners, lon_corners; n_lon, n_lat)

Compute exact planar overlap weights between one quadrilateral footprint and
regular grid cells in fractional grid-index coordinates. The returned matrix has
shape `(n_lon, n_lat)` and sums to one when the footprint is fully inside the
domain.
"""
function exact_footprint_weights(lat_corners, lon_corners; n_lon::Integer, n_lat::Integer)
    poly = _footprint_polygon(lat_corners, lon_corners)
    total_area = polygon_area(poly)
    weights = zeros(Float32, n_lon, n_lat)
    total_area <= eps(Float64) && return weights

    min_col = max(1, floor(Int, minimum(first, poly)))
    max_col = min(n_lon, floor(Int, maximum(first, poly)))
    min_row = max(1, floor(Int, minimum(last, poly)))
    max_row = min(n_lat, floor(Int, maximum(last, poly)))

    for col in min_col:max_col, row in min_row:max_row
        clipped = _clip_polygon_to_cell(poly, col, row)
        weights[col, row] = Float32(polygon_area(clipped) / total_area)
    end

    weights
end

"""
    footprint_weights(method, lat_corners, lon_corners; n_lon, n_lat)

Return per-cell footprint weights for a single footprint using the requested
geometry method. This is a lightweight geometry-level API used by tests and by
future exact-intersection gridders.
"""
function footprint_weights(method::ExactIntersectionGridding, lat_corners, lon_corners;
                           n_lon::Integer, n_lat::Integer)
    exact_footprint_weights(lat_corners, lon_corners; n_lon=n_lon, n_lat=n_lat)
end

function footprint_weights(method::SubpixelGridding, lat_corners, lon_corners;
                           n_lon::Integer, n_lat::Integer)
    n = method.n_oversample
    n === nothing && error("Subpixel footprint_weights requires an explicit n_oversample")

    points = zeros(Float32, n, n, 2)
    ix = zeros(Int32, n^2)
    iy = zeros(Int32, n^2)
    lats0 = zeros(Float32, n)
    lons0 = zeros(Float32, n)
    lats1 = zeros(Float32, n)
    lons1 = zeros(Float32, n)

    lat_sorted, lon_sorted = sort_corners_ccw(Float32.(lat_corners),
                                              Float32.(lon_corners))
    compute_subpixels!(points, lat_sorted, lon_sorted, n,
                       lats0, lons0, lats1, lons1)
    floor_indices!(ix, iy, points)

    weights = zeros(Float32, n_lon, n_lat)
    fac = Float32(1 / n^2)
    @inbounds for i in eachindex(ix)
        row = ix[i]
        col = iy[i]
        if 1 <= col <= n_lon && 1 <= row <= n_lat
            weights[col, row] += fac
        end
    end
    weights
end

function footprint_weights(::CenterPointGridding, lat_corners, lon_corners;
                           n_lon::Integer, n_lat::Integer)
    center_lat = sum(lat_corners) / length(lat_corners)
    center_lon = sum(lon_corners) / length(lon_corners)
    row = floor(Int, center_lat)
    col = floor(Int, center_lon)

    weights = zeros(Float32, n_lon, n_lat)
    if 1 <= col <= n_lon && 1 <= row <= n_lat
        weights[col, row] = 1.0f0
    end
    weights
end
