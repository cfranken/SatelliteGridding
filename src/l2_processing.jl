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

# Spherical grids locate each subpoint independently. Their subdivision
# trait selects the geometry and sampling buffers.
function _process_l2_file!(filepath::String, config::DataSourceConfig,
                           spec::AbstractSphericalGridSpec, grid_data, grid_std, grid_weights,
                           compute_std, n_override, method::AbstractGriddingMethod,
                           nc_vars, fill_attrib, progress; metadata=nothing)
    _with_l2_observations(filepath, config, spec, method, nc_vars, fill_attrib,
                          progress; metadata=metadata) do observations
        _accumulate_l2_spherical!(_spherical_subdivision(spec), grid_data, grid_std,
                                  grid_weights, compute_std, spec, observations, n_override)
    end
    nothing
end

function _accumulate_l2_spherical!(::Val{:lonlat}, data, m2, weights, compute_std,
                                   spec, observations, n_override)
    n = n_override === nothing ? 4 : n_override
    points = zeros(_spherical_subpixel_type(spec), n, n, 2)
    lats0, lons0, lats1, lons1 = zeros(n), zeros(n), zeros(n), zeros(n)
    values = observations.values
    n_pixels, n_vars = size(values)
    _accumulate_spherical_footprint!(data, m2, weights, compute_std, spec,
                                     observations.lat_bounds, observations.lon_bounds,
                                     values, n_pixels, n_vars, n,
                                     points, lats0, lons0, lats1, lons1)
end

function _accumulate_l2_spherical!(::Val{:unit_vectors}, data, m2, weights, compute_std,
                                   spec, observations, n_override)
    capacity = n_override === nothing ? _HEALPIX_MAX_OVERSAMPLE : n_override
    cells = Vector{Tuple{Int,Int}}(undef, capacity^2)
    _accumulate_unit_vector_footprints!(data, m2, weights, compute_std, spec,
                                       observations.lat_bounds, observations.lon_bounds,
                                       observations.values, n_override, cells)
end

function _process_l2_file_ka!(filepath::String, config::DataSourceConfig,
                              spec::GridSpec{T}, backend, grid_sum, grid_weights,
                              n_override, method::AbstractGriddingMethod,
                              nc_vars, fill_attrib, progress;
                              metadata=nothing, spherical_ka=nothing, reference=nothing,
                              pass=nothing) where {T}
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

function _process_l2_file_ka!(filepath::String, config::DataSourceConfig,
                              spec::AbstractSphericalGridSpec, backend, grid_sum, grid_weights,
                              n_override, method::SubpixelGridding,
                              nc_vars, fill_attrib, progress;
                              spherical_ka, reference, metadata=nothing, pass=nothing)
    # Screen native Float32/Float64 science values before narrowing for a
    # Float32 accumulator. Finite values outside its range must reach the
    # kernel so their cells get fill for the non-finite moments.
    _with_l2_observations(filepath, config, spec, method, nc_vars, fill_attrib,
                          progress; metadata, value_type=Float64) do observations
        _accumulate_spherical_batch!(backend, grid_sum, grid_weights, pass,
            spherical_ka.locator, spherical_ka.geometry_precision, observations.lat_bounds,
            observations.lon_bounds, observations.values, n_override, reference)
    end
    return nothing
end

"""
    _accumulate_l2_std_pass!(files, config, spec, backend, m2, residual, weights,
                            mean, n_override, method, spherical_ka;
                            reference, bounds, progress, keep_shifted=false)

Replay the first-pass files around each cell mean, then correct the mean
and central M2 using the accumulated linear and squared deviations. A replay
error invalidates the window because that file already contributed to the mean.
With `keep_shifted=true`, return the reference-shifted mean in `mean` and its
low-order correction in `residual` for a later rolling-window merge.
"""
function _accumulate_l2_std_pass!(files, config, spec, backend, m2, residual, weights,
                                  mean, n_override, method, spherical_ka;
                                  reference, bounds, progress, keep_shifted=false)
    _finalize_spherical_mean!(backend, mean, weights, reference, bounds)
    fill!(m2, zero(eltype(m2)))
    fill!(residual, zero(eltype(residual)))
    p_files = _progress_bar(length(files), "  STD files: ", progress)
    for filepath in files
        try
            _process_l2_file_ka!(filepath, config, spec, backend, m2, weights,
                                 n_override, method, nothing, false, p_files;
                                 spherical_ka, reference, pass=(mean, residual))
        catch e
            error("Error replaying file $filepath for standard deviation: $(sprint(showerror, e))")
        end
    end
    _correct_spherical_moments_ka!(backend, _SPHERICAL_WORKGROUP_SIZE)(
        mean, m2, residual, weights, reference.device, Val(keep_shifted); ndrange=length(weights))
    KernelAbstractions.synchronize(backend)
    return nothing
end
