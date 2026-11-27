# Read and prepare observations independently of the accumulation backend.
_has_l2_bounds(config::DataSourceConfig) =
    haskey(config.basic, "lat_bnd") && haskey(config.basic, "lon_bnd")

_has_l2_centers(config::DataSourceConfig) =
    haskey(config.basic, "lat") && haskey(config.basic, "lon")

_has_l2_radius(config::DataSourceConfig) =
    haskey(config.basic, "radius") || haskey(config.options, "radius")

_validate(spec, method::AbstractGriddingMethod) = nothing

function _validate(spec::AbstractSphericalGridSpec, ::CircularFootprintGridding)
    error("CircularFootprintGridding is unsupported for $(typeof(spec)). Use SubpixelGridding (quad footprints).")
end

_use_radius_l2_geometry(config, ::AbstractGriddingMethod) = false
_use_radius_l2_geometry(config, ::CircularFootprintGridding) =
    _has_l2_radius(config) && !_has_l2_bounds(config)

function _validate_l2_geometry_config(config::DataSourceConfig, method::AbstractGriddingMethod)
    has_lat = haskey(config.basic, "lat")
    has_lon = haskey(config.basic, "lon")
    has_lat == has_lon ||
        error("grid_l2 requires both basic.lat and basic.lon, or neither so centers can be derived from bounds")
    _validate_l2_coordinates(config, method)
    if !has_lat
        @warn "basic.lat/basic.lon not configured; deriving footprint centers from corner coordinates"
    end
    nothing
end

function _validate_l2_coordinates(config, ::SubpixelGridding)
    _has_l2_bounds(config) ||
        error("grid_l2 requires basic.lat_bnd and basic.lon_bnd for quadrilateral footprint gridding")
    nothing
end

function _validate_l2_coordinates(config, ::CircularFootprintGridding)
    _has_l2_bounds(config) && return nothing
    _has_l2_radius(config) ||
        error("CircularFootprintGridding requires either basic.lat_bnd/basic.lon_bnd or center basic.lat/basic.lon plus basic.radius or [circle] radius")
    _has_l2_centers(config) ||
        error("Radius-based CircularFootprintGridding requires basic.lat and basic.lon center variables")
    nothing
end

function _read_or_derive_l2_centers(fin, config::DataSourceConfig,
                                    lat_bnd::AbstractMatrix,
                                    lon_bnd::AbstractMatrix)
    _has_l2_centers(config) && return _read_l2_centers(fin, config)

    vec(mean(lat_bnd, dims=2)), vec(mean(lon_bnd, dims=2))
end

function _read_l2_centers(fin, config::DataSourceConfig)
    _has_l2_centers(config) ||
        error("Radius-based circular footprint gridding requires basic.lat and basic.lon")
    vec(_read_l2_geolocation(fin, config.basic["lat"])),
    vec(_read_l2_geolocation(fin, config.basic["lon"]))
end

function _read_l2_radius(fin, config::DataSourceConfig, n::Integer)
    if haskey(config.basic, "radius")
        radius = vec(read_nc_variable(fin, config.basic["radius"]))
        length(radius) == n ||
            error("Radius variable '$(config.basic["radius"])' has $(length(radius)) values, expected $n")
        return radius
    end

    haskey(config.options, "radius") ||
        error("Radius-based circular footprint gridding requires basic.radius or [circle] radius")
    fill(Float64(config.options["radius"]), n)
end

function _radius_to_degrees(radius::Real, lat::Real, unit::AbstractString)
    unit_lc = lowercase(strip(unit))
    if unit_lc in ("degree", "degrees", "deg")
        r = Float64(radius)
    elseif unit_lc in ("kilometer", "kilometers", "km")
        r = Float64(radius) / 111.32
    elseif unit_lc in ("meter", "meters", "m")
        r = Float64(radius) / 111_320.0
    else
        error("Unknown [circle] radius_unit '$unit'. Use degrees, km, or m.")
    end

    lon_scale = max(abs(cosd(Float64(lat))), 1e-6)
    r, r / lon_scale
end

function _circle_bounds_from_radius(lat_center, lon_center, radius,
                                    config::DataSourceConfig)
    unit = string(get(config.options, "radius_unit", "degrees"))
    n = length(lat_center)
    lat_bnd = zeros(Float64, n, 4)
    lon_bnd = zeros(Float64, n, 4)

    @inbounds for i in 1:n
        if !isfinite(radius[i]) || radius[i] < 0
            lat_bnd[i, :] .= NaN
            lon_bnd[i, :] .= NaN
            continue
        end
        rlat, rlon = _radius_to_degrees(radius[i], lat_center[i], unit)
        lat_bnd[i, 1] = lat_center[i] - rlat
        lat_bnd[i, 2] = lat_center[i] - rlat
        lat_bnd[i, 3] = lat_center[i] + rlat
        lat_bnd[i, 4] = lat_center[i] + rlat
        lon_bnd[i, 1] = lon_center[i] - rlon
        lon_bnd[i, 2] = lon_center[i] + rlon
        lon_bnd[i, 3] = lon_center[i] + rlon
        lon_bnd[i, 4] = lon_center[i] - rlon
    end

    lat_bnd, lon_bnd
end

function _read_l2_geometry(fin, config::DataSourceConfig,
                           footprint_method::AbstractGriddingMethod, spec::AbstractGridSpec)
    if _use_radius_l2_geometry(config, footprint_method)
        lat_center, lon_center = _read_l2_centers(fin, config)
        radius = _read_l2_radius(fin, config, length(lat_center))
        lat_bnd, lon_bnd = _circle_bounds_from_radius(lat_center, lon_center,
                                                      radius, config)
        return lat_center, lon_center, lat_bnd, lon_bnd
    end

    lat_bnd = _read_l2_geolocation(fin, config.basic["lat_bnd"]; bounds=true)
    lon_bnd = _read_l2_geolocation(fin, config.basic["lon_bnd"]; bounds=true)

    _order_l2_corners!(spec, lat_bnd, lon_bnd)
    lat_center, lon_center = _read_or_derive_l2_centers(fin, config,
                                                        lat_bnd, lon_bnd)
    lat_center, lon_center, lat_bnd, lon_bnd
end

_order_l2_corners!(::RectangularGridSpec, lat, lon) = sort_corners_ccw!(lat, lon)
_order_l2_corners!(::AbstractSphericalGridSpec, lat, lon) = nothing

const _POLAR_LATITUDE_TOLERANCE = 1e-4 # degrees
const _SPHERICAL_LONGITUDE_LIMIT = 720 # degrees

@inline _clamp_polar_latitude(lat) = lat
@inline function _clamp_polar_latitude(lat::T) where {T<:AbstractFloat}
    return 90 < abs(lat) <= 90 + _POLAR_LATITUDE_TOLERANCE ? copysign(T(90), lat) : lat
end

_validate_l2_geolocation!(valid, ::RectangularGridSpec, lat, lon, lat_c, lon_c) = nothing

"""
Clamp small polar overshoots and reject spherical soundings whose centre or
corners fall outside the latitude and longitude domains in physical degrees.
`valid` already excludes nonfinite coordinates.
"""
function _validate_l2_geolocation!(valid, ::AbstractSphericalGridSpec, lat, lon, lat_c, lon_c)
    @inbounds for i in eachindex(valid)
        valid[i] || continue
        lat[i] = _clamp_polar_latitude(lat[i])
        valid[i] = (-90 <= lat[i] <= 90) & (-_SPHERICAL_LONGITUDE_LIMIT <= lon[i] <= _SPHERICAL_LONGITUDE_LIMIT)
        for j in axes(lat_c, 2)
            lat_c[i, j] = _clamp_polar_latitude(lat_c[i, j])
            valid[i] &= -90 <= lat_c[i, j] <= 90
            valid[i] &= -_SPHERICAL_LONGITUDE_LIMIT <= lon_c[i, j] <= _SPHERICAL_LONGITUDE_LIMIT
        end
    end
    nothing
end

# Use the same input preparation for rectangular and spherical grids and KA backends.
function _with_l2_observations(f, filepath, config, spec::AbstractGridSpec{T}, method,
                               nc_vars, fill_attrib, progress;
                               metadata=nothing, value_type=T) where {T}
    Dataset(filepath) do fin
        _capture_input_metadata!(metadata, fin)
        if fill_attrib && nc_vars !== nothing
            attributes = get(get(config.options, "metadata", Dict()),
                             "variable_attributes", DEFAULT_VARIABLE_ATTRIBUTES)
            copy_variable_attributes!(nc_vars, fin, config.grid_vars; attributes)
        end
        observations = _prepare_l2_observations(fin, config, spec, method, progress, filepath; value_type)
        observations === nothing || f(observations)
    end
    nothing
end

_l2_file_in_bounds(::AbstractSphericalGridSpec, lat, lon) = true
_l2_file_in_bounds(spec::RectangularGridSpec, lat, lon) =
    any(i -> in_bounds(spec, lat[i], lon[i]), eachindex(lat, lon))

function _prepare_l2_observations(fin, config, spec::AbstractGridSpec{T}, method,
                                  progress, filepath; value_type=T) where {T}
    lat_center, lon_center, lat_bounds, lon_bounds = _read_l2_geometry(fin, config, method, spec)
    if !_l2_file_in_bounds(spec, lat_center, lon_center)
        ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, 0)])
        return nothing
    end
    idx = apply_filters(fin, config, lat_bounds, lon_bounds, spec)
    ProgressMeter.next!(progress; showvalues=[(:File, filepath), (:N_pixels, length(idx))])
    isempty(idx) && return nothing

    # Materialize only selected rows in the shared matrix, rather than an
    # all-soundings × all-variables matrix followed by another full copy.
    values = Matrix{value_type}(undef, length(idx), length(config.grid_vars))
    for (column, path) in enumerate(Base.values(config.grid_vars))
        data = _read_l2_values(fin, path, value_type)
        length(data) == size(lat_bounds, 1) ||
            throw(DimensionMismatch("Variable '$path' has $(length(data)) values; expected $(size(lat_bounds, 1)) soundings"))
        copyto!(view(values, :, column), view(data, idx))
    end
    lat = lat_center[idx]
    lon = lon_center[idx]
    lat_c = lat_bounds[idx, :]
    lon_c = lon_bounds[idx, :]
    finite = vec(all(isfinite, values; dims=2)) .&
             vec(all(isfinite, lat_c; dims=2)) .& vec(all(isfinite, lon_c; dims=2)) .&
             isfinite.(lat) .& isfinite.(lon)
    _validate_l2_geolocation!(finite, spec, lat, lon, lat_c, lon_c)
    any(finite) || return nothing
    if !all(finite)
        values = values[finite, :]
        lat, lon = lat[finite], lon[finite]
        lat_c, lon_c = lat_c[finite, :], lon_c[finite, :]
    end
    (; values, lat, lon, lat_bounds=lat_c, lon_bounds=lon_c)
end

function _rectangular_indices(observations, spec::RectangularGridSpec, method)
    n_lat, n_lon = length(spec.lat), length(spec.lon)
    lat_scale = spec.lat_max - spec.lat_min
    lon_scale = spec.lon_max - spec.lon_min
    lat_index(lat) = ((lat .- spec.lat_min) / lat_scale * n_lat) .+ 1
    lon_index(lon) = ((lon .- spec.lon_min) / lon_scale * n_lon) .+ 1
    ilat = lat_index(observations.lat_bounds)
    ilon = lon_index(observations.lon_bounds)
    centers = _rectangular_centers(method, observations, lat_index, lon_index)
    ilat, ilon, centers
end

_rectangular_centers(::SubpixelGridding, obs, lat_index, lon_index) = nothing
_rectangular_centers(::CircularFootprintGridding, obs, lat_index, lon_index) =
    (lat_index(obs.lat), lon_index(obs.lon))

function _l2_oversampling(ilat, ilon, override)
    override === nothing || return override
    lat_extents = maximum(ilat; dims=2) .- minimum(ilat; dims=2)
    lon_extents = maximum(ilon; dims=2) .- minimum(ilon; dims=2)
    compute_n_oversample(median(max.(lat_extents, lon_extents)), 1.0)
end
