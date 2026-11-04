# Backend-specific accumulation of prepared L2 observations.
function _process_l2_file!(filepath::String, config::DataSourceConfig,
                           spec::GridSpec{T}, grid_data, grid_std, grid_weights,
                           compute_std, n_override, method::AbstractGriddingMethod,
                           nc_vars, fill_attrib, progress; metadata=nothing) where {T}
    _with_l2_observations(filepath, config, spec, method, nc_vars, fill_attrib,
                          progress; metadata=metadata) do observations
        ilat, ilon, ilat_center, ilon_center = _rectangular_indices(observations, spec)
        n = _l2_oversampling(ilat, ilon, n_override)
        values = observations.values
        n_pixels, n_vars = size(values)
        if method isa CircularFootprintGridding
            accumulate_circular_footprint!(grid_data, grid_std, grid_weights,
                                           compute_std, ilat_center, ilon_center,
                                           ilat, ilon, values, n_pixels, n_vars, n)
        else
            points = zeros(T, n, n, 2)
            ix, iy = zeros(Int32, n^2), zeros(Int32, n^2)
            lats0, lons0, lats1, lons1 = zeros(n), zeros(n), zeros(n), zeros(n)
            accumulate_footprint!(grid_data, grid_std, grid_weights, compute_std,
                                  ilat, ilon, values, n_pixels, n_vars, n,
                                  points, ix, iy, lats0, lons0, lats1, lons1)
        end
    end
    nothing
end

# Cubed-sphere accumulation subdivides in lon/lat space before projecting each
# point; rectangular affine indices and their extent heuristic do not apply.
function _process_l2_file!(filepath::String, config::DataSourceConfig,
                           spec::CubedSphereGridSpec{T}, grid_data, grid_std, grid_weights,
                           compute_std, n_override, method::AbstractGriddingMethod,
                           nc_vars, fill_attrib, progress; metadata=nothing) where {T}
    _with_l2_observations(filepath, config, spec, method, nc_vars, fill_attrib,
                          progress; metadata=metadata) do observations
        n = n_override === nothing ? 4 : n_override
        points = zeros(T, n, n, 2)
        lats0, lons0, lats1, lons1 = zeros(n), zeros(n), zeros(n), zeros(n)
        values = observations.values
        n_pixels, n_vars = size(values)
        accumulate_footprint_cs!(grid_data, grid_std, grid_weights, compute_std,
                                 spec, observations.lat_bounds, observations.lon_bounds,
                                 values, n_pixels, n_vars, n,
                                 points, lats0, lons0, lats1, lons1)
    end
    nothing
end

function _process_l2_file_ka!(filepath::String, config::DataSourceConfig,
                              spec::GridSpec{T}, backend, grid_sum, grid_weights,
                              n_override, method::AbstractGriddingMethod,
                              nc_vars, fill_attrib, progress; metadata=nothing) where {T}
    _with_l2_observations(filepath, config, spec, method, nc_vars, fill_attrib,
                          progress; metadata=metadata) do observations
        ilat, ilon, ilat_center, ilon_center =
            map(indices -> T.(indices), _rectangular_indices(observations, spec))
        n = _l2_oversampling(ilat, ilon, n_override)
        values = observations.values
        n_vars = size(values, 2)
        if method isa CircularFootprintGridding
            accumulate_circular_batch!(backend, grid_sum, grid_weights,
                                       ilat_center, ilon_center, ilat, ilon,
                                       values, n, n_vars)
        else
            accumulate_batch!(backend, grid_sum, grid_weights, ilat, ilon,
                              values, n, n_vars)
        end
    end
    nothing
end
