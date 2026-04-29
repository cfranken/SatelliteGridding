"""
    apply_filters(dataset, config::DataSourceConfig, lat_bounds, lon_bounds,
                  grid_spec::GridSpec) -> Vector{Int}

Apply spatial bounding box and quality filters, returning indices of pixels that pass
all criteria.

The spatial filter requires all 4 corners to be within the grid bounds and the footprint
to not cross the dateline (lon span < 50°). Quality filters from `config.filters` are
applied as additional criteria.
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
        data = read_nc_variable(dataset, rule.variable)
        if rule.op === :lt
            bool_add .+= (data .< rule.lo)
        elseif rule.op === :gt
            bool_add .+= (data .> rule.lo)
        elseif rule.op === :eq
            bool_add .+= (data .== rule.lo)
        elseif rule.op === :between
            bool_add .+= ((data .> rule.lo) .& (data .< rule.hi))
        end
        n_criteria += 1
    end

    findall(bool_add .== n_criteria)
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
