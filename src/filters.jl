@inline function _passes_quality_filter(value, rule::FilterRule)
    (ismissing(value) || !isfinite(value)) && return false
    rule.op === :lt && return value < rule.lo
    rule.op === :gt && return value > rule.lo
    rule.op === :eq && return value == rule.lo
    rule.lo < value < rule.hi
end

function _spatial_mask(spec::RectangularGridSpec, lat_bounds, lon_bounds)
    min_lat = minimum(lat_bounds, dims=2)
    max_lat = maximum(lat_bounds, dims=2)
    min_lon = minimum(lon_bounds, dims=2)
    max_lon = maximum(lon_bounds, dims=2)

    # Add the criteria so missing coordinates propagate even when another
    # spatial or quality criterion fails.
    ((min_lat[:, 1] .> spec.lat_min) .+
     (max_lat[:, 1] .< spec.lat_max) .+
     (min_lon[:, 1] .> spec.lon_min) .+
     (max_lon[:, 1] .< spec.lon_max) .+
     ((max_lon[:, 1] .- min_lon[:, 1]) .< 50)) .== 5
end

_spatial_mask(::AbstractSphericalGridSpec, lat_bounds, lon_bounds) = trues(size(lat_bounds, 1))

_combine_quality(spatial::Bool, quality::Bool) = spatial & quality
_combine_quality(::Missing, ::Bool) = missing

"""
    apply_filters(dataset, config, lat_bounds, lon_bounds, grid_spec) -> Vector{Int}

Return the indices of pixels passing all configured quality filters. Rectangular
grids also require all four footprint corners strictly inside the domain and a
longitude span below 50°. Spherical grids cover the entire sphere, including
dateline and polar footprints. Input preparation checks their coordinate validity
before accumulation.
"""
function apply_filters(dataset, config::DataSourceConfig,
                       lat_bounds::AbstractMatrix, lon_bounds::AbstractMatrix,
                       spec::AbstractGridSpec)::Vector{Int}
    mask = _spatial_mask(spec, lat_bounds, lon_bounds)
    for rule in config.filters
        data = _read_l2_values(dataset, rule.variable, Float64)
        mask .= _combine_quality.(mask, _passes_quality_filter.(data, Ref(rule)))
    end
    findall(mask)
end

"""
    apply_center_filters(lat, lon, grid_spec::RectangularGridSpec) -> Vector{CartesianIndex}

Apply spatial bounding box filter for center-coordinate gridding (no footprint bounds).
Returns CartesianIndices of pixels within the grid bounds.
"""
function apply_center_filters(lat::AbstractArray, lon::AbstractArray,
                              grid_spec::RectangularGridSpec)
    findall((lat .> grid_spec.lat_min) .& (lat .< grid_spec.lat_max) .&
            (lon .> grid_spec.lon_min) .& (lon .< grid_spec.lon_max))
end
