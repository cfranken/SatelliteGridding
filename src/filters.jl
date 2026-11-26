@inline function _passes_quality_filter(value, rule::FilterRule)
    (ismissing(value) || !isfinite(value)) && return false
    rule.op === :lt && return value < rule.lo
    rule.op === :gt && return value > rule.lo
    rule.op === :eq && return value == rule.lo
    rule.lo < value < rule.hi
end

"""
    apply_filters(dataset, config::DataSourceConfig, lat_bounds, lon_bounds,
                  grid_spec::GridSpec) -> Vector{Int}

Apply spatial bounding box and quality filters, returning indices of pixels that pass
all criteria.

The spatial filter requires all four footprint corners to lie strictly inside
`(lat_min, lat_max) × (lon_min, lon_max)` and the footprint not to cross the dateline
(`lon span < 50°`); footprints with any corner on or outside the domain boundary are
excluded before accumulation. Quality filters from `config.filters` are applied as
additional criteria.
"""
function apply_filters(dataset, config::DataSourceConfig,
                       lat_bounds::AbstractMatrix, lon_bounds::AbstractMatrix,
                       grid_spec::GridSpec)::Vector{Int}
    min_lat = minimum(lat_bounds, dims=2)
    max_lat = maximum(lat_bounds, dims=2)
    min_lon = minimum(lon_bounds, dims=2)
    max_lon = maximum(lon_bounds, dims=2)

    # Spatial bounding box + dateline check
    bool_add = ((min_lat[:, 1] .> grid_spec.lat_min) .+
                (max_lat[:, 1] .< grid_spec.lat_max) .+
                (min_lon[:, 1] .> grid_spec.lon_min) .+
                (max_lon[:, 1] .< grid_spec.lon_max) .+
                ((max_lon[:, 1] .- min_lon[:, 1]) .< 50))
    n_criteria = 5

    for rule in config.filters
        data = _read_l2_values(dataset, rule.variable, Float64)
        bool_add .+= _passes_quality_filter.(data, Ref(rule))
        n_criteria += 1
    end

    findall(bool_add .== n_criteria)
end

"""
    apply_filters(dataset, config, lat_bounds, lon_bounds,
                  grid_spec::AbstractSphericalGridSpec) -> Vector{Int}

Cubed-sphere and HEALPix grids cover the whole sphere, so this method applies
only the quality filters from `config.filters`. A valid high-latitude or
near-pole footprint can span a large longitude range, including across the
dateline, so there is no spatial bounding box or dateline span check.
Input preparation checks that the centre and all corners are valid before
the sounding can contribute to the accumulator.
"""
function apply_filters(dataset, config::DataSourceConfig,
                       lat_bounds::AbstractMatrix, lon_bounds::AbstractMatrix,
                       grid_spec::AbstractSphericalGridSpec)::Vector{Int}
    n = size(lat_bounds, 1)
    bool_add = zeros(Int, n)
    n_criteria = 0

    for rule in config.filters
        data = _read_l2_values(dataset, rule.variable, Float64)
        bool_add .+= _passes_quality_filter.(data, Ref(rule))
        n_criteria += 1
    end

    n_criteria == 0 ? collect(1:n) : findall(bool_add .== n_criteria)
end

"""
    apply_center_filters(lat, lon, grid_spec::GridSpec) -> Vector{CartesianIndex}

Apply spatial bounding box filter for center-coordinate gridding (no footprint bounds).
Returns CartesianIndices of pixels within the grid bounds.
"""
function apply_center_filters(lat::AbstractArray, lon::AbstractArray,
                              grid_spec::GridSpec)
    findall((lat .> grid_spec.lat_min) .& (lat .< grid_spec.lat_max) .&
            (lon .> grid_spec.lon_min) .& (lon .< grid_spec.lon_max))
end
