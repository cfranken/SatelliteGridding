# HEALPix base-face offsets in the polar/equatorial projection of Górski et al.
# Faces 0–3, 4–7, and 8–11 are the north, equatorial, and south belts.
const _HEALPIX_RING_OFFSETS = (2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
const _HEALPIX_LONGITUDE_OFFSETS = (1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7)
const _HEALPIX_FACE_COUNT = length(_HEALPIX_RING_OFFSETS)
const _HEALPIX_TRANSITION_Z = 2 / 3
const _HEALPIX_POLAR_COLATITUDE = 0.01  # healpy's angular projection crossover, in radians
# The legacy angular projection loses precision near its polar crossover.
const _HEALPIX_TIE_TOLERANCE = 128eps(Float64)
# Largest power of two for which the cell count fits in Int64.
const _HEALPIX_MAX_NSIDE = Int64(1) << 29
const _HEALPIX_EARTH_RADIUS = 6.371e6  # metres
const _HEALPIX_COORD_CHUNK = 1 << 16  # cells per coordinate write
const _HEALPIX_VERTICES = 4
const _HEALPIX_SAMPLES_PER_CELL = 4
const _HEALPIX_MIN_OVERSAMPLE = 4
const _HEALPIX_MAX_OVERSAMPLE = 24

function _healpix_nside(nside)
    nside isa Integer && !(nside isa Bool) && ispow2(nside) &&
        1 <= nside <= _HEALPIX_MAX_NSIDE ||
        throw(ArgumentError("HEALPix nside must be a power of two between 1 and $_HEALPIX_MAX_NSIDE"))
    return Int64(nside)
end

"""
    nside2npix(nside) -> Int64

Return the number of HEALPix cells, `12*nside^2`, for a NESTED grid.
`nside` must be a power of two between 1 and 2^29.
"""
function nside2npix(nside::Integer)
    side = _healpix_nside(nside)
    return _HEALPIX_FACE_COUNT * side^2
end

"""
    healpix_parent(pix) -> Int64

Return the zero-based NESTED index one level above `pix`. The given pixel
must belong to a level above the base resolution (`nside > 1`).
"""
healpix_parent(pix::Integer) = Int64(pix) >> 2

# Masks that spread or collect bits between alternating positions of a
# 64-bit Morton code. NESTED stores x in the even bits and y in the odd bits.
const _HEALPIX_MORTON_MASKS = (0x5555555555555555, 0x3333333333333333,
    0x0f0f0f0f0f0f0f0f, 0x00ff00ff00ff00ff, 0x0000ffff0000ffff, 0x00000000ffffffff)

@inline function _healpix_spread_bits(x::Int64)
    m1, m2, m4, m8, m16, _ = _HEALPIX_MORTON_MASKS
    bits = reinterpret(UInt64, x)
    bits = (bits | (bits << 16)) & m16
    bits = (bits | (bits << 8)) & m8
    bits = (bits | (bits << 4)) & m4
    bits = (bits | (bits << 2)) & m2
    return (bits | (bits << 1)) & m1
end

@inline function _healpix_compact_bits(bits::UInt64)
    m1, m2, m4, m8, m16, m32 = _HEALPIX_MORTON_MASKS
    bits &= m1
    bits = (bits | (bits >> 1)) & m2
    bits = (bits | (bits >> 2)) & m4
    bits = (bits | (bits >> 4)) & m8
    bits = (bits | (bits >> 8)) & m16
    return Int64((bits | (bits >> 16)) & m32)
end

"""
    ang2pix_nest(nside, lon_deg, lat_deg) -> Int64

Locate a point on a HEALPix grid using the zero-based NESTED convention of
Górski et al. (2005, ApJ 622, 759). Longitude is periodic, with edge ties
evaluated on the canonical branch `[0, 360)` degrees. Latitude is spherical,
in degrees within `[-90, 90]`.
"""
@inline function ang2pix_nest(nside::Integer, lon_deg::Real, lat_deg::Real)
    side = _healpix_nside(nside)
    return _ang2pix_nest(side, trailing_zeros(side), Float64(lon_deg), Float64(lat_deg))
end

@inline function _ang2pix_nest(nside::Int64, level::Int, lon::Float64, lat::Float64)
    isfinite(lon) && -90 <= lat <= 90 ||
        throw(ArgumentError("HEALPix requires finite longitude and latitude in [-90, 90]"))
    return _vec2pix_nest(nside, level, _healpix_angular_projection(nside, lon, lat)...)
end

@inline function _healpix_angular_projection(nside, lon, lat)
    # Preserve healpy's operation order at exact angular pixel edges, including
    # the longitude-dependent face at an angular pole.
    theta = (pi / 2) - deg2rad(lat)
    z = cos(theta)
    za = abs(z)
    quadrant = deg2rad(_canonical_longitude(lon)) * (2 / pi)
    quadrant = quadrant == 4.0 ? 0.0 : quadrant
    scale = if za <= _HEALPIX_TRANSITION_Z
        0.0
    elseif theta < _HEALPIX_POLAR_COLATITUDE || theta > pi - _HEALPIX_POLAR_COLATITUDE
        nside * sin(theta) / _spherical_sqrt((1 + za) / 3)
    else
        nside * _spherical_sqrt(3 * (1 - za))
    end
    return z, quadrant, scale
end

"""
    vec2pix_nest(nside, x, y, z) -> Int64

Locate a finite unit vector in zero-based HEALPix NESTED order.
The components share a precision of Float32 or Float64. At an exact pole
with zero transverse components the vector convention uses longitude zero.
Use Float64 for very fine grids. Above roughly `nside = 2^21`, Float32
pixel errors can exceed one cell.
"""
@inline function vec2pix_nest(nside::Integer, x::F, y::F, z::F) where {F<:Union{Float32,Float64}}
    side = _healpix_nside(nside)
    return _vec2pix_nest(side, trailing_zeros(side), x, y, z, Val(false))
end

@inline function _vec2pix_nest(nside::Int64, level::Int, x::Float64, y::Float64,
                              z::Float64, ::Val{P}) where {P}
    q = (iszero(x) && iszero(y)) ? 0.0 : atan(y, x) * (2 / pi)
    quadrant = q < 0.0 ? q + 4.0 : q
    quadrant = quadrant == 4.0 ? 0.0 : quadrant
    za = abs(z)
    scale = za <= _HEALPIX_TRANSITION_Z ? 0.0 :
        nside * _spherical_sqrt(3.0 * ((x * x + y * y) / (1.0 + za)))
    if P && _healpix_near_edge(nside, z, quadrant, scale)
        # The sequential footprint path historically round-tripped through
        # angular coordinates. Retain its decisions on numerically tied edges.
        horizontal = hypot(x, y)
        lon = iszero(horizontal) ? 0.0 : atand(y, x)
        lat = atand(z, horizontal)
        return _vec2pix_nest(nside, level, _healpix_angular_projection(nside, lon, lat)...)
    end
    return _vec2pix_nest(nside, level, z, quadrant, scale)
end

@inline function _healpix_quadrant_fraction(x::Float32, y::Float32)
    # Rotate to the first quadrant using only coordinate swaps and sign flips.
    if y > 0 && x <= 0
        quadrant, u, v = 1, y, -x
    elseif x < 0 && y <= 0
        quadrant, u, v = 2, -x, -y
    elseif y < 0
        quadrant, u, v = 3, -y, x
    else
        quadrant, u, v = 0, x, y
    end
    fraction = iszero(v) ? 0f0 : min(prevfloat(1f0), atan(v, u) * Float32(2 / pi))
    return quadrant, fraction
end

@inline function _vec2pix_nest(nside::Int64, level::Int, x::Float32, y::Float32,
                              z::Float32, ::Val)
    quadrant, fraction = _healpix_quadrant_fraction(x, y)
    za = abs(z)
    equatorial = za <= Float32(_HEALPIX_TRANSITION_Z)
    if equatorial
        # Add whole quadrants and nside/2 in integer space so their magnitude
        # cannot round away the low bits of the fractional longitude.
        offset = (quadrant << level) + (nside >> 1)
        along_lon = Float32(nside) * fraction + Float32(nside & 1) * 0.5f0
        along_lat = Float32(nside) * (z * 0.75f0)
        ascending = offset + Base.unsafe_trunc(Int64, floor(along_lon - along_lat))
        descending = offset + Base.unsafe_trunc(Int64, floor(along_lon + along_lat))
    else
        # Resolve polar rings even if Float32 rounds abs(z) to one.
        scale = Float32(nside) * _spherical_sqrt(3f0 * ((x * x + y * y) / (1f0 + za)))
        ascending = Base.unsafe_trunc(Int64, fraction * scale)
        descending = Base.unsafe_trunc(Int64, (1f0 - fraction) * scale)
    end
    return _healpix_nest_from_edges(nside, level, z, equatorial, quadrant, ascending, descending)
end

@inline function _healpix_edges(nside, z::F, quadrant::F, scale::F) where {F}
    if abs(z) <= F(_HEALPIX_TRANSITION_Z)
        along_lon = F(nside) * (F(0.5) + quadrant)
        along_lat = F(nside) * (z * F(0.75))
        return along_lon - along_lat, along_lon + along_lat
    end
    fraction = quadrant - F(min(3, Base.unsafe_trunc(Int64, quadrant)))
    return fraction * scale, (one(F) - fraction) * scale
end

@inline function _healpix_near_edge(nside, z, quadrant, scale)
    tolerance = _HEALPIX_TIE_TOLERANCE * nside
    abs(abs(z) - _HEALPIX_TRANSITION_Z) <= _HEALPIX_TIE_TOLERANCE && return true
    a, b = _healpix_edges(nside, z, quadrant, scale)
    return abs(a - round(a)) <= tolerance || abs(b - round(b)) <= tolerance
end

@inline function _vec2pix_nest(nside::Int64, level::Int, z::F, quadrant::F, scale::F) where {F}
    a, b = _healpix_edges(nside, z, quadrant, scale)
    ascending = Base.unsafe_trunc(Int64, a)
    descending = Base.unsafe_trunc(Int64, b)
    face = min(3, Base.unsafe_trunc(Int64, quadrant))
    return _healpix_nest_from_edges(nside, level, z, abs(z) <= F(_HEALPIX_TRANSITION_Z),
                                  face, ascending, descending)
end

@inline function _healpix_nest_from_edges(nside, level, z, equatorial, face, ascending, descending)
    last = nside - 1
    if equatorial
        # The ascending and descending edge lines form a square grid.
        asc_face = ascending >> level
        desc_face = descending >> level
        # At the longitude wrap, rounding can put an edge on sector 4.
        # Wrap the sector within its belt before adding the belt offset.
        face = if asc_face == desc_face
            4 + (asc_face & 3)
        elseif asc_face < desc_face
            asc_face & 3
        else
            8 + (desc_face & 3)
        end
        x = descending & last
        y = last - (ascending & last)
    else
        ascending = min(last, ascending)
        descending = min(last, descending)
        if z >= 0
            x, y = last - descending, last - ascending
        else
            x, y = ascending, descending
            face += 8
        end
    end
    local_id = _healpix_spread_bits(x) | (_healpix_spread_bits(y) << 1)
    return (face << (2 * level)) + reinterpret(Int64, local_id)
end

@inline function _healpix_face_xy(level::Int, pix::Int64)
    face = pix >> (2 * level)
    local_id = UInt64(pix - (face << (2 * level)))
    return face + 1, _healpix_compact_bits(local_id), _healpix_compact_bits(local_id >> 1)
end

# Convert fractional face coordinates (0–1 along each face edge) to degrees.
# nr is the latitude ring's radius in the planar HEALPix projection.
@inline function _healpix_face_lonlat(face, x, y)
    ring = _HEALPIX_RING_OFFSETS[face] - x - y
    if ring < 1 || ring > 3
        north = ring < 1
        nr = north ? ring : 4 - ring
        delta = nr^2 / 3
        # Use the transverse component to retain precision next to the poles.
        lat = atand(sqrt(delta * (2 - delta)), 1 - delta)
        lat = north ? 90 - lat : lat - 90
        nr == 0 && return 0.0, lat
    else
        nr = 1.0
        lat = asind((2 - ring) * _HEALPIX_TRANSITION_Z)
    end
    lon = 45 * (_HEALPIX_LONGITUDE_OFFSETS[face] + (x - y) / nr)
    return mod(lon, 360.0), lat
end

function _check_healpix_pixel(nside::Int64, pix::Integer)
    0 <= pix < _HEALPIX_FACE_COUNT * nside^2 ||
        throw(ArgumentError("HEALPix pixel index is out of range for nside=$nside"))
    return Int64(pix)
end

"""
    pix2ang_nest(nside, pix) -> (lon_deg, lat_deg)

Return the spherical coordinates in degrees of the centre of zero-based
NESTED pixel `pix`. Longitude is in `[0, 360)`.
"""
function pix2ang_nest(nside::Integer, pix::Integer)
    side = _healpix_nside(nside)
    return _pix2ang_nest(side, trailing_zeros(side), _check_healpix_pixel(side, pix))
end

@inline function _pix2ang_nest(nside::Int64, level::Int, pix::Int64)
    face, x, y = _healpix_face_xy(level, pix)
    return _healpix_face_lonlat(face, (x + 0.5) / nside, (y + 0.5) / nside)
end

"""
    healpix_corners(nside, pix) -> (lons, lats)

Return the four vertex longitudes and latitudes in degrees for zero-based
NESTED pixel `pix`, as tuples in north, west, south, east order, matching
`healpy.boundaries(step=1)`. Longitude is in `[0, 360)`.
"""
function healpix_corners(nside::Integer, pix::Integer)
    side = _healpix_nside(nside)
    return _healpix_corners(side, trailing_zeros(side), _check_healpix_pixel(side, pix))
end

@inline function _healpix_corners(nside::Int64, level::Int, pix::Int64)
    face, x, y = _healpix_face_xy(level, pix)
    x0, x1 = x / nside, (x + 1) / nside
    y0, y1 = y / nside, (y + 1) / nside
    lon_n, lat_n = _healpix_face_lonlat(face, x1, y1)
    lon_w, lat_w = _healpix_face_lonlat(face, x0, y1)
    lon_s, lat_s = _healpix_face_lonlat(face, x0, y0)
    lon_e, lat_e = _healpix_face_lonlat(face, x1, y0)
    return (lon_n, lon_w, lon_s, lon_e), (lat_n, lat_w, lat_s, lat_e)
end

"""
    HEALPixGridSpec(; nside, T=Float32, cell_bounds=false)

Whole-sphere HEALPix grid in NESTED order with equal solid angle per cell.
The resolution is
`nside = 2^level`. Accumulators fold the cells into `(nside^2, 12)` so linear
index `pix + 1` addresses NESTED cell `pix`. Sequential accumulators and
finalized slices use `T`. Select KA geometry and accumulator precision with
`grid_l2`. Float32 accumulators (used on Metal) can report standard deviations
below roughly one ulp of the cell mean (`σ/|mean| ≲ 1e-7`) as zero. Use
Float64 accumulators (the CUDA and KA CPU default) for small relative variability.
Mean-only KA runs skip refinement. In either accumulator precision `A`, values
with magnitude below about `eps(A) * abs(reference - value)` can round away
when shifted by the per-variable reference (`eps(Float64)` ≈ `2e-16`,
`eps(Float32)` ≈ `1e-7`). Enable `compute_std` to recover those values.
Sequential projections and output coordinates use Float64.
Cell centres and optional vertex bounds are computed when writing NetCDF.

The grid uses input geodetic latitude directly as spherical latitude,
consistent with the cubed sphere. True WGS84 cell areas vary by about
1.35% from equator to pole (−0.45% at the equator and +0.90% at the poles
relative to the nominal spherical `cell_area`).
"""
struct HEALPixGridSpec{T<:AbstractFloat} <: AbstractSphericalGridSpec{T}
    nside::Int64
    level::Int
    cell_bounds::Bool

    function HEALPixGridSpec(; nside, T::Type{<:AbstractFloat}=Float32, cell_bounds::Bool=false)
        side = _healpix_nside(nside)
        return new{T}(side, trailing_zeros(side), cell_bounds)
    end
end

grid_shape(spec::HEALPixGridSpec) = (spec.nside^2, _HEALPIX_FACE_COUNT)
payload_dims(::HEALPixGridSpec) = ("cell",)
_spherical_subdivision(::HEALPixGridSpec) = Val(:unit_vectors)

# Angular form of sqrt(4πR²/npix). The common radius also scales the
# angular footprint span, so it cancels in the oversampling ratio.
_healpix_cell_spacing(spec::HEALPixGridSpec) = _device_cell_spacing(_device_locator(spec), Float64)

struct _HEALPixLocator
    nside::Int64
    level::Int
end

_device_locator(spec::HEALPixGridSpec) = _HEALPixLocator(spec.nside, spec.level)

@inline function locate(locator::_HEALPixLocator, x, y, z)
    return 1 + _vec2pix_nest(locator.nside, locator.level, x, y, z, Val(true))
end

@inline function _device_cell_spacing(locator::_HEALPixLocator, ::Type{F}) where {F}
    return sqrt(F(4pi) / F(_HEALPIX_FACE_COUNT * locator.nside^2))
end

@inline function _unit_vector_samples(corners, cell_spacing::F) where {F}
    max_chord_squared = zero(F)
    @inbounds for i in 1:_HEALPIX_VERTICES, j in (i + 1):_HEALPIX_VERTICES
        a, b = corners[i], corners[j]
        chord_squared = abs2(a[1] - b[1]) + abs2(a[2] - b[2]) + abs2(a[3] - b[3])
        max_chord_squared = max(max_chord_squared, chord_squared)
    end
    span = F(2) * asin(min(one(F), _spherical_sqrt(max_chord_squared) / F(2)))
    samples = Base.unsafe_trunc(Int, ceil(F(_HEALPIX_SAMPLES_PER_CELL) * span / cell_spacing))
    return clamp(samples, _HEALPIX_MIN_OVERSAMPLE, _HEALPIX_MAX_OVERSAMPLE)
end

@inline in_bounds(::HEALPixGridSpec, lat::Real, lon::Real) =
    -90 <= lat <= 90 && isfinite(lon)

@inline function _spherical_cell(spec::HEALPixGridSpec, lat::Real, lon::Real)
    pix = _ang2pix_nest(spec.nside, spec.level, Float64(lon), Float64(lat))
    face = pix >> (2 * spec.level)
    return pix - (face << (2 * spec.level)) + 1, face + 1
end

function write_coordinates!(ds, spec::HEALPixGridSpec)
    npix = nside2npix(spec.nside)
    defDim(ds, "cell", npix)
    cell_dims = payload_dims(spec)
    chunk = min(npix, _HEALPIX_COORD_CHUNK)
    cell_chunks = [chunk]
    ids = defVar(ds, "cell_ids", Int64, cell_dims; chunksizes=cell_chunks,
        attrib=["standard_name" => "healpix_index", "units" => "1",
                "grid_name" => "healpix", "level" => spec.level,
                "indexing_scheme" => "nested", "ellipsoid" => "sphere"])
    lon_attrs = ["standard_name" => "longitude", "units" => "degrees_east"]
    lat_attrs = ["standard_name" => "latitude", "units" => "degrees_north"]
    if spec.cell_bounds
        defDim(ds, "nv", _HEALPIX_VERTICES)
        push!(lon_attrs, "bounds" => "lon_bounds")
        push!(lat_attrs, "bounds" => "lat_bounds")
        bounds_dims = ("nv", cell_dims...)
        bounds_chunk = [_HEALPIX_VERTICES, chunk]
        bounds_attrs = ["comment" => "Corner vertices for plotting. HEALPix edges are curved; great-circle polygons differ in area by about 10% for cells touching the poles."]
        lon_bounds = defVar(ds, "lon_bounds", Float64, bounds_dims; chunksizes=bounds_chunk, attrib=bounds_attrs)
        lat_bounds = defVar(ds, "lat_bounds", Float64, bounds_dims; chunksizes=bounds_chunk, attrib=bounds_attrs)
    end
    lon = defVar(ds, "lon", Float64, cell_dims; chunksizes=cell_chunks, attrib=lon_attrs)
    lat = defVar(ds, "lat", Float64, cell_dims; chunksizes=cell_chunks, attrib=lat_attrs)
    mapping = defVar(ds, "healpix", Int32, ();
        attrib=["grid_mapping_name" => "healpix", "refinement_level" => spec.level,
                "indexing_scheme" => "nested",
                "earth_radius" => _HEALPIX_EARTH_RADIUS])
    area = defVar(ds, "cell_area", Float64, ();
        attrib=["standard_name" => "cell_area", "units" => "m2",
                "comment" => "Cells subtend equal solid angle on the sphere. Geodetic latitude is used as spherical latitude, as for the cubed sphere. True WGS84 areas vary by about 1.35% from equator to pole (−0.45% at the equator, +0.90% at the poles relative to cell_area)."])
    mapping[] = Int32(0)
    area[] = 4pi * _HEALPIX_EARTH_RADIUS^2 / npix

    lons, lats = Vector{Float64}(undef, chunk), Vector{Float64}(undef, chunk)
    if spec.cell_bounds
        corner_lons = Matrix{Float64}(undef, _HEALPIX_VERTICES, chunk)
        corner_lats = similar(corner_lons)
    end
    for first in 1:chunk:npix
        cells = first:min(first + chunk - 1, npix)
        columns = 1:length(cells)
        @inbounds for (j, cell) in enumerate(cells)
            lons[j], lats[j] = _pix2ang_nest(spec.nside, spec.level, cell - 1)
        end
        ids[cells] = (first - 1):(last(cells) - 1)
        lon[cells] = @view lons[columns]
        lat[cells] = @view lats[columns]
        if spec.cell_bounds
            @inbounds for (j, cell) in enumerate(cells)
                xs, ys = _healpix_corners(spec.nside, spec.level, cell - 1)
                for vertex in 1:_HEALPIX_VERTICES
                    corner_lons[vertex, j] = xs[vertex]
                    corner_lats[vertex, j] = ys[vertex]
                end
            end
            lon_bounds[:, cells] = @view corner_lons[:, columns]
            lat_bounds[:, cells] = @view corner_lats[:, columns]
        end
    end

    ds.attrib["Conventions"] = "CF-1.13"
    ds.attrib["grid"] = "healpix"
    ds.attrib["nside"] = spec.nside
    ds.attrib["healpix_level"] = spec.level
    ds.attrib["healpix_ordering"] = "nested"
    return nothing
end
