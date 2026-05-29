# ---------------------------------------------------------------------------
# Cubed-sphere (GEOS/GCHP) output grid.
#
# The closed-form gnomonic projection and its inverse `lonlat_to_panel_xy`
# below are PORTED VERBATIM (math unchanged) from the AtmosTransportModel
# package:
#
#   /home/cfranken/code/gitHub/AtmosTransportModel/src/Grids/CubedSphereMesh.jl
#   commit 6522d85 ("Add closed-form cubed-sphere inverse projection
#                    lonlat_to_panel_xy")
#
# Only the self-contained projection chain is reproduced here — the mesh,
# panel connectivity, spherical-area metrics, and tangent-basis helpers from
# the source are intentionally omitted because the gridder needs only the
# forward map (to build cell-center / corner lon-lat arrays for the NetCDF
# coordinates) and the inverse map (to scatter observations onto cells).
#
# "Cubed sphere" is not a single geometry: GEOS-IT/GEOS-FP (and therefore
# GCHP when driven by GMAO met) use the GMAO equal-distance gnomonic grid with
# GEOS-native panel ordering and a -10 deg longitude shift; older synthetic
# targets used an equiangular gnomonic construction. Both are supported via the
# coordinate-law / center-law / panel-convention contract.
#
# References:
#   Putman & Lin (2007) — FV3 cubed-sphere grid
#   Martin et al. (2022, GMD) — GCHP v13
# ---------------------------------------------------------------------------

abstract type AbstractCubedSphereCoordinateLaw end
abstract type AbstractCubedSphereCenterLaw end
abstract type AbstractCubedSpherePanelConvention end
abstract type AbstractCubedSphereDefinition end

"""
    EquiangularGnomonic <: AbstractCubedSphereCoordinateLaw

Classical equiangular gnomonic face coordinate law. For edge index
`s = 1, ..., Nc + 1` the local angle is `α_s = -π/4 + (s - 1) π / (2Nc)` and the
tangent-plane coordinate is `ξ_s = tan(α_s)`. Useful for synthetic targets; it
is **not** the native GEOS-IT/GEOS-FP coordinate law.
"""
struct EquiangularGnomonic <: AbstractCubedSphereCoordinateLaw end

raw"""
    GMAOEqualDistanceGnomonic <: AbstractCubedSphereCoordinateLaw

GMAO/GEOS equal-distance gnomonic edge coordinate law (`grid_type = 0` in the
GEOS/FV and ESMF code path) used by GEOS-IT C180 and GEOS-FP C720. With
`r = 1/√3` and `α₀ = sin⁻¹(r)`, for edge index `s = 1, ..., Nc + 1`:

```math
β_s = -\frac{α_0}{Nc}(Nc + 2 - 2s), \qquad b_s = \tan(β_s)\cos(α_0),
```

and the tangent-plane coordinate is `ξ_s = b_s / r`. This non-uniform edge law
produces the characteristic GEOS C180 meridional spacing (~0.42° at panel edges,
~0.55° near panel centers).
"""
struct GMAOEqualDistanceGnomonic <: AbstractCubedSphereCoordinateLaw end

"""
    AngularMidpointCenter <: AbstractCubedSphereCenterLaw

Cell-center law that evaluates the continuous face map at the logical cell
midpoint `x_c(i,j) = x(i + 1/2, j + 1/2)`. Matches the historical synthetic /
equiangular behavior; not how GEOS writes native center coordinates.
"""
struct AngularMidpointCenter <: AbstractCubedSphereCenterLaw end

raw"""
    FourCornerNormalizedCenter <: AbstractCubedSphereCenterLaw

GEOS/FV/ESMF `cell_center2` center law: the cell center is the normalized sum of
the four unit corner vectors,
`v_c = (v₁ + v₂ + v₃ + v₄) / ‖v₁ + v₂ + v₃ + v₄‖`. GEOS-IT/GEOS-FP `lons`/`lats`
arrays are written from this law.
"""
struct FourCornerNormalizedCenter <: AbstractCubedSphereCenterLaw end

"""
    GnomonicPanelConvention <: AbstractCubedSpherePanelConvention

Classical gnomonic panel numbering: 1-4 equatorial, 5 north pole, 6 south pole.
"""
struct GnomonicPanelConvention <: AbstractCubedSpherePanelConvention end

"""
    GEOSNativePanelConvention <: AbstractCubedSpherePanelConvention

Panel numbering/orientation used by native GEOS-FP / GEOS-IT cubed-sphere files:
panels 1-2 equatorial, 3 north pole, 4-5 equatorial, 6 south pole. Describes
panel storage/order only; the GMAO `-10°` shift and equal-distance law live in
[`GMAOCubedSphereDefinition`](@ref).
"""
struct GEOSNativePanelConvention <: AbstractCubedSpherePanelConvention end

"""
    CubedSphereDefinition(coordinate_law, center_law, panel_convention;
                          longitude_offset_deg=0, tag=:custom)

Complete horizontal cubed-sphere geometry contract: how logical face edges are
placed on the gnomonic cube (`coordinate_law`), how cell centers are derived
(`center_law`), how the six faces are ordered/oriented (`panel_convention`), and
a final rigid z-axis rotation (`longitude_offset_deg`). For native GMAO files use
[`GMAOCubedSphereDefinition`](@ref).
"""
struct CubedSphereDefinition{L<:AbstractCubedSphereCoordinateLaw,
                             C<:AbstractCubedSphereCenterLaw,
                             P<:AbstractCubedSpherePanelConvention} <: AbstractCubedSphereDefinition
    coordinate_law::L
    center_law::C
    panel_convention::P
    longitude_offset_deg::Float64
    tag::Symbol
end

function CubedSphereDefinition(law::L, center::C, convention::P;
                               longitude_offset_deg::Real=0,
                               tag::Symbol=:custom) where {
                               L<:AbstractCubedSphereCoordinateLaw,
                               C<:AbstractCubedSphereCenterLaw,
                               P<:AbstractCubedSpherePanelConvention}
    return CubedSphereDefinition{L,C,P}(law, center, convention,
                                        Float64(longitude_offset_deg), tag)
end

"""
    EquiangularCubedSphereDefinition(; convention=GnomonicPanelConvention(),
                                       longitude_offset_deg=0)

Legacy synthetic/equiangular definition: `EquiangularGnomonic` corners plus
`AngularMidpointCenter` centers.
"""
EquiangularCubedSphereDefinition(;
    convention::AbstractCubedSpherePanelConvention=GnomonicPanelConvention(),
    longitude_offset_deg::Real=0) =
    CubedSphereDefinition(EquiangularGnomonic(), AngularMidpointCenter(), convention;
                          longitude_offset_deg=longitude_offset_deg,
                          tag=:equiangular_gnomonic)

"""
    GMAOCubedSphereDefinition(; convention=GEOSNativePanelConvention(),
                                longitude_offset_deg=-10)

Native GMAO/GEOS cubed-sphere definition used by GEOS-IT C180 / GEOS-FP C720 (and
GCHP when driven by GMAO met): `GMAOEqualDistanceGnomonic` edges,
`FourCornerNormalizedCenter` centers, `GEOSNativePanelConvention`, and a final
`-10°` longitude shift. This is the default for [`CubedSphereGridSpec`](@ref).
"""
GMAOCubedSphereDefinition(;
    convention::AbstractCubedSpherePanelConvention=GEOSNativePanelConvention(),
    longitude_offset_deg::Real=-10) =
    CubedSphereDefinition(GMAOEqualDistanceGnomonic(), FourCornerNormalizedCenter(),
                          convention; longitude_offset_deg=longitude_offset_deg,
                          tag=:gmao_equal_distance)

coordinate_law(def::CubedSphereDefinition) = def.coordinate_law
center_law(def::CubedSphereDefinition) = def.center_law
panel_convention(def::CubedSphereDefinition) = def.panel_convention
longitude_offset_deg(def::CubedSphereDefinition) = def.longitude_offset_deg
cs_definition_tag(def::CubedSphereDefinition) = def.tag

coordinate_law_tag(::EquiangularGnomonic) = "equiangular_gnomonic"
coordinate_law_tag(::GMAOEqualDistanceGnomonic) = "gmao_equal_distance_gnomonic"
center_law_tag(::AngularMidpointCenter) = "angular_midpoint"
center_law_tag(::FourCornerNormalizedCenter) = "four_corner_normalized"

# ---------------------------------------------------------------------------
# Gnomonic projection (forward)
# ---------------------------------------------------------------------------

"""
    _gnomonic_xyz(ξ, η, panel) -> (x, y, z)

Map local tangent-plane coordinates `(ξ, η)` to unit-sphere Cartesian
`(x, y, z)` via the gnomonic (central) projection for `panel`. The
normalisation `d = 1/√(1 + ξ² + η²)` projects onto the unit sphere. Panels 1-4
are the equatorial belt; 5, 6 are the polar caps (see source for the per-panel
axis table).
"""
@inline function _gnomonic_xyz(ξ::FT, η::FT, panel::Int) where {FT}
    d = one(FT) / sqrt(one(FT) + ξ^2 + η^2)   # gnomonic normalisation
    if panel == 1
        return (d, ξ * d, η * d)              # +x face
    elseif panel == 2
        return (-ξ * d, d, η * d)             # +y face
    elseif panel == 3
        return (-d, -ξ * d, η * d)            # −x face
    elseif panel == 4
        return (ξ * d, -d, η * d)             # −y face
    elseif panel == 5
        return (-η * d, ξ * d, d)             # +z (north pole)
    else
        return (η * d, ξ * d, -d)            # −z (south pole)
    end
end

"Convert unit-sphere Cartesian to `(lon, lat)` in degrees, lon in [0, 360)."
@inline function _xyz_to_lonlat(x, y, z)
    lon = atand(y, x)
    lat = asind(z / sqrt(x^2 + y^2 + z^2))
    lon < 0 && (lon += 360)
    return lon, lat
end

@inline function _normalize3(x, y, z)
    n = sqrt(x^2 + y^2 + z^2)
    invn = inv(max(n, eps(typeof(n))))
    return (x * invn, y * invn, z * invn)
end

@inline function _rotate_z_lon_offset(x::FT, y::FT, z::FT, offset_deg::Real) where {FT}
    θ = FT(deg2rad(offset_deg))
    c = cos(θ)
    s = sin(θ)
    return (c * x - s * y, s * x + c * y, z)
end

@inline function _edge_tangent_coordinate(::EquiangularGnomonic, s::Real, Nc::Int,
                                          ::Type{FT}) where {FT}
    α = -FT(π) / 4 + (FT(s) - one(FT)) * (FT(π) / (2 * FT(Nc)))
    return tan(α)
end

@inline function _edge_tangent_coordinate(::GMAOEqualDistanceGnomonic, s::Real, Nc::Int,
                                          ::Type{FT}) where {FT}
    r = inv(sqrt(FT(3)))
    α0 = asin(r)
    β = -(α0 / FT(Nc)) * (FT(Nc) + FT(2) - FT(2) * FT(s))
    b = tan(β) * cos(α0)
    return b / r
end

@inline function _panel_xyz(::GnomonicPanelConvention, ξ::FT, η::FT, panel::Int) where {FT}
    return _gnomonic_xyz(ξ, η, panel)
end

# GEOS-FP/GEOS-IT native panel ordering 1,2,north,4,5,south. Panels 4 and 5 are
# stored 90° CW rotated relative to the gnomonic order, and panel 3 (north) with
# a quarter-turn. Validated against GEOSIT.20211202.A3dyn.C180.nc (see source).
@inline function _panel_xyz(::GEOSNativePanelConvention, ξ::FT, η::FT, panel::Int) where {FT}
    ξg, ηg, gpanel = if panel == 1
        (ξ, η, 1)
    elseif panel == 2
        (ξ, η, 2)
    elseif panel == 3
        (-η, ξ, 5)
    elseif panel == 4
        (η, -ξ, 3)
    elseif panel == 5
        (η, -ξ, 4)
    elseif panel == 6
        (ξ, η, 6)
    else
        throw(ArgumentError("invalid GEOS native panel id $panel"))
    end
    return _gnomonic_xyz(ξg, ηg, gpanel)
end

@inline function _panel_xyz(def::CubedSphereDefinition, ξ::FT, η::FT, panel::Int) where {FT}
    x, y, z = _panel_xyz(panel_convention(def), ξ, η, panel)
    offset = longitude_offset_deg(def)
    return iszero(offset) ? (x, y, z) : _rotate_z_lon_offset(x, y, z, offset)
end

@inline function _continuous_panel_xyz(def::CubedSphereDefinition, Nc::Int,
                                       s::Real, t::Real, panel::Int, ::Type{FT}) where {FT}
    law = coordinate_law(def)
    ξ = _edge_tangent_coordinate(law, s, Nc, FT)
    η = _edge_tangent_coordinate(law, t, Nc, FT)
    return _panel_xyz(def, ξ, η, panel)
end

@inline function _corner_xyz(def::CubedSphereDefinition, Nc::Int,
                             i::Integer, j::Integer, panel::Int, ::Type{FT}) where {FT}
    return _continuous_panel_xyz(def, Nc, i, j, panel, FT)
end

@inline function _cell_center_xyz(def::CubedSphereDefinition, ::AngularMidpointCenter,
                                  Nc::Int, i::Integer, j::Integer, panel::Int,
                                  ::Type{FT}) where {FT}
    return _continuous_panel_xyz(def, Nc, FT(i) + FT(0.5), FT(j) + FT(0.5), panel, FT)
end

@inline function _cell_center_xyz(def::CubedSphereDefinition, ::FourCornerNormalizedCenter,
                                  Nc::Int, i::Integer, j::Integer, panel::Int,
                                  ::Type{FT}) where {FT}
    v1 = _corner_xyz(def, Nc, i, j, panel, FT)
    v2 = _corner_xyz(def, Nc, i + 1, j, panel, FT)
    v3 = _corner_xyz(def, Nc, i + 1, j + 1, panel, FT)
    v4 = _corner_xyz(def, Nc, i, j + 1, panel, FT)
    return _normalize3(v1[1] + v2[1] + v3[1] + v4[1],
                       v1[2] + v2[2] + v3[2] + v4[2],
                       v1[3] + v2[3] + v3[3] + v4[3])
end

@inline function _cell_center_xyz(def::CubedSphereDefinition, Nc::Int,
                                  i::Integer, j::Integer, panel::Int, ::Type{FT}) where {FT}
    return _cell_center_xyz(def, center_law(def), Nc, i, j, panel, FT)
end

# ---------------------------------------------------------------------------
# Inverse projection: (lon, lat) → (panel, s_frac, t_frac)
# ---------------------------------------------------------------------------

"Inverse of [`_gnomonic_xyz`](@ref): unit `(x,y,z)` → `(ξg, ηg, gpanel)`."
@inline function _gnomonic_inverse(x::FT, y::FT, z::FT) where {FT}
    ax = abs(x)
    ay = abs(y)
    az = abs(z)
    if ax >= ay && ax >= az
        if x > 0
            return (y / x, z / x, 1)
        else
            return (y / x, -z / x, 3)
        end
    elseif ay >= az
        if y > 0
            return (-x / y, z / y, 2)
        else
            return (-x / y, -z / y, 4)
        end
    else
        if z > 0
            return (y / z, -x / z, 5)
        else
            return (-y / z, -x / z, 6)
        end
    end
end

"Undo the panel re-mapping `_panel_xyz(convention, …)` applied before the gnomonic map."
@inline function _undo_panel_convention(::GnomonicPanelConvention,
                                        ξg::FT, ηg::FT, gpanel::Int) where {FT}
    return (ξg, ηg, gpanel)
end

@inline function _undo_panel_convention(::GEOSNativePanelConvention,
                                        ξg::FT, ηg::FT, gpanel::Int) where {FT}
    # panel 3 → gpanel 5 with (ξg,ηg)=(-η,ξ)  ⇒ inverse (ξ,η)=(ηg,-ξg)
    # panel 4 → gpanel 3 with (ξg,ηg)=( η,-ξ) ⇒ inverse (ξ,η)=(-ηg,ξg)
    # panel 5 → gpanel 4 with (ξg,ηg)=( η,-ξ) ⇒ inverse (ξ,η)=(-ηg,ξg)
    if gpanel == 1
        return (ξg, ηg, 1)
    elseif gpanel == 2
        return (ξg, ηg, 2)
    elseif gpanel == 5
        return (ηg, -ξg, 3)
    elseif gpanel == 3
        return (-ηg, ξg, 4)
    elseif gpanel == 4
        return (-ηg, ξg, 5)
    else
        return (ξg, ηg, 6)
    end
end

@inline function _inverse_edge_tangent_coordinate(::EquiangularGnomonic, ξ::Real, Nc::Int,
                                                  ::Type{FT}) where {FT}
    α = atan(FT(ξ))
    return one(FT) + (α + FT(π) / 4) * (FT(2 * Nc) / FT(π))
end

@inline function _inverse_edge_tangent_coordinate(::GMAOEqualDistanceGnomonic, ξ::Real, Nc::Int,
                                                  ::Type{FT}) where {FT}
    r = inv(sqrt(FT(3)))
    α0 = asin(r)
    b = FT(ξ) * r
    β = atan(b / cos(α0))
    return FT(Nc + 2) / FT(2) + β * FT(Nc) / (FT(2) * α0)
end

"""
    lonlat_to_panel_xy(def::CubedSphereDefinition, Nc, lon_deg, lat_deg, FT)
        -> (panel::Int, s_frac, t_frac)

Closed-form inverse of the cubed-sphere forward map for a single `(lon, lat)` in
degrees. Returns the panel id `panel ∈ 1..6` and continuous edge indices
`s_frac`, `t_frac` on the panel's `[1, Nc+1)²` parameter square (half-open at the
upper edge). The integer cell index is `floor(Int, s_frac)` (x) and
`floor(Int, t_frac)` (y), guaranteed in `1..Nc` for every finite `(lon, lat)` —
upper-edge points are nudged to `prevfloat(Nc+1)` and sub-`1` round-off is raised
to `1`, the standard half-open-cell convention. The whole sphere is covered.

Honors the definition's panel convention, coordinate law, and longitude offset.
"""
@inline function lonlat_to_panel_xy(def::CubedSphereDefinition, Nc::Int,
                                    lon_deg::Real, lat_deg::Real, ::Type{FT}) where {FT}
    coslat = cosd(FT(lat_deg))
    x = coslat * cosd(FT(lon_deg))
    y = coslat * sind(FT(lon_deg))
    z = sind(FT(lat_deg))

    offset = longitude_offset_deg(def)
    if !iszero(offset)
        x, y, z = _rotate_z_lon_offset(x, y, z, -offset)
    end

    ξg, ηg, gpanel = _gnomonic_inverse(x, y, z)
    ξ, η, panel = _undo_panel_convention(panel_convention(def), ξg, ηg, gpanel)

    law = coordinate_law(def)
    s = _inverse_edge_tangent_coordinate(law, ξ, Nc, FT)
    t = _inverse_edge_tangent_coordinate(law, η, Nc, FT)

    s_max = prevfloat(FT(Nc + 1))
    t_max = prevfloat(FT(Nc + 1))
    s = clamp(s, one(FT), s_max)
    t = clamp(t, one(FT), t_max)
    return (panel, s, t)
end

# ---------------------------------------------------------------------------
# CubedSphereGridSpec — AbstractGridSpec implementation
# ---------------------------------------------------------------------------

"""
    CubedSphereGridSpec{T,D} <: AbstractGridSpec{T}

Cubed-sphere output grid with `Nc` cells per panel edge and 6 panels (`6·Nc²`
cells total). The default `definition` is [`GMAOCubedSphereDefinition`](@ref) —
the GEOS/GCHP convention (equal-distance gnomonic, GEOS-native panels, `-10°`
longitude shift), e.g. `CubedSphereGridSpec(Nc=360)` for C360.

# Accumulator folding
Internally the cube `(Nc, Nc, 6)` is *folded* into a 2D accumulator of shape
`grid_shape(spec) = (Nc, 6·Nc)` so the existing rectangular gridder machinery
(`_welford_update!`, `finalize_mean!`, `accumulate_center!`) is reused unchanged:
cell `(i, j, panel)` maps to `col = i`, `row = (panel-1)·Nc + j`.
[`to_fractional_index`](@ref) returns the folded `(col, row)`; the NetCDF writer
un-folds via `reshape(slice, Nc, Nc, 6)` back to `(Xdim, Ydim, nf)`.

# Fields
- `Nc::Int`: cells per panel edge
- `definition::D`: cubed-sphere geometry contract
- `centers_lon`, `centers_lat`: `(Nc, Nc, 6)` cell-center coords (degrees, lon in [0,360))
- `corners_lon`, `corners_lat`: `(Nc+1, Nc+1, 6)` cell-corner coords (degrees)
"""
struct CubedSphereGridSpec{T<:AbstractFloat,D<:CubedSphereDefinition} <: AbstractGridSpec{T}
    Nc::Int
    definition::D
    centers_lon::Array{T,3}
    centers_lat::Array{T,3}
    corners_lon::Array{T,3}
    corners_lat::Array{T,3}
end

"""
    CubedSphereGridSpec(; Nc, definition=GMAOCubedSphereDefinition(), T=Float32)

Construct a cubed-sphere grid spec, precomputing the cell-center and cell-corner
lon/lat arrays from `definition`. Use `Nc=360` for a GCHP-style C360 grid.
"""
function CubedSphereGridSpec(; Nc::Integer,
                              definition::CubedSphereDefinition=GMAOCubedSphereDefinition(),
                              T::Type{<:AbstractFloat}=Float32)
    Nc > 0 || throw(ArgumentError("Nc must be positive, got $Nc"))
    Nci = Int(Nc)

    centers_lon = Array{T,3}(undef, Nci, Nci, 6)
    centers_lat = Array{T,3}(undef, Nci, Nci, 6)
    @inbounds for p in 1:6, j in 1:Nci, i in 1:Nci
        x, y, z = _cell_center_xyz(definition, Nci, i, j, p, T)
        lon, lat = _xyz_to_lonlat(x, y, z)
        centers_lon[i, j, p] = T(lon)
        centers_lat[i, j, p] = T(lat)
    end

    corners_lon = Array{T,3}(undef, Nci + 1, Nci + 1, 6)
    corners_lat = Array{T,3}(undef, Nci + 1, Nci + 1, 6)
    @inbounds for p in 1:6, j in 1:(Nci + 1), i in 1:(Nci + 1)
        x, y, z = _corner_xyz(definition, Nci, i, j, p, T)
        lon, lat = _xyz_to_lonlat(x, y, z)
        corners_lon[i, j, p] = T(lon)
        corners_lat[i, j, p] = T(lat)
    end

    return CubedSphereGridSpec{T,typeof(definition)}(Nci, definition,
                                                     centers_lon, centers_lat,
                                                     corners_lon, corners_lat)
end

# ---------------------------------------------------------------------------
# AbstractGridSpec interface — CubedSphereGridSpec implementations
# ---------------------------------------------------------------------------

"""
    to_fractional_index(spec::CubedSphereGridSpec, lat, lon) -> (col, row)

Project `(lat, lon)` (degrees) onto continuous indices in the **folded**
`(Nc, 6·Nc)` accumulator layout: `col = s_frac ∈ [1, Nc+1)` and
`row = (panel-1)·Nc + t_frac`, so `floor(Int, col) = i ∈ 1..Nc` and
`floor(Int, row) = (panel-1)·Nc + j ∈ 1..6·Nc`.
"""
@inline function to_fractional_index(spec::CubedSphereGridSpec{T}, lat::Real, lon::Real) where {T}
    panel, s, t = lonlat_to_panel_xy(spec.definition, spec.Nc, lon, lat, T)
    col = s
    row = T((panel - 1) * spec.Nc) + t
    return (col, row)
end

# Whole-sphere grid: every finite (lat, lon) lands in exactly one cell.
@inline in_bounds(::CubedSphereGridSpec, lat::Real, lon::Real) =
    isfinite(lat) && isfinite(lon)

grid_shape(spec::CubedSphereGridSpec) = (spec.Nc, 6 * spec.Nc)

payload_dims(::CubedSphereGridSpec) = ("Xdim", "Ydim", "nf")

function write_coordinates!(ds, spec::CubedSphereGridSpec)
    Nc = spec.Nc
    # Dimension names and ordering follow native GEOS/GCHP cubed-sphere output so
    # the files open directly in Panoply, xarray, and GCHP tooling. NetCDF stores
    # dimensions in the reverse of how they are declared here (column-major
    # Julia), so declaring (Xdim, Ydim, nf) yields CDL `(nf, Ydim, Xdim)` —
    # exactly the GCHP layout. Corners use the (Nc+1) `XCdim`/`YCdim` dims.
    defDim(ds, "Xdim", Nc)
    defDim(ds, "Ydim", Nc)
    defDim(ds, "nf", 6)
    defDim(ds, "XCdim", Nc + 1)
    defDim(ds, "YCdim", Nc + 1)

    # 1D "fake" GrADS-compatibility axes in degrees (GCHP convention): the
    # panel-1 cell-center longitudes/latitudes along each logical axis. The real
    # georeferencing is the 2D `lons`/`lats` below (`coordinates = "lons lats"`).
    ds_x = defVar(ds, "Xdim", Float64, ("Xdim",),
                  attrib=["long_name" => "Fake Longitude for GrADS Compatibility",
                          "units" => "degrees_east"])
    ds_y = defVar(ds, "Ydim", Float64, ("Ydim",),
                  attrib=["long_name" => "Fake Latitude for GrADS Compatibility",
                          "units" => "degrees_north"])
    ds_nf = defVar(ds, "nf", Int32, ("nf",),
                   attrib=["long_name" => "cubed-sphere face", "axis" => "e"])
    ds_x[:] = Float64.(@view spec.centers_lon[:, 1, 1])
    ds_y[:] = Float64.(@view spec.centers_lat[1, :, 1])
    ds_nf[:] = Int32.(1:6)

    # 2D cell-center aux coordinates (degrees). Order (Xdim, Ydim, nf) → CDL
    # (nf, Ydim, Xdim), matching GCHP `lons`/`lats`.
    ds_lons = defVar(ds, "lons", Float64, ("Xdim", "Ydim", "nf"),
                     attrib=["units" => "degrees_east", "long_name" => "longitude"])
    ds_lats = defVar(ds, "lats", Float64, ("Xdim", "Ydim", "nf"),
                     attrib=["units" => "degrees_north", "long_name" => "latitude"])
    ds_lons[:, :, :] = Float64.(spec.centers_lon)
    ds_lats[:, :, :] = Float64.(spec.centers_lat)

    ds_clons = defVar(ds, "corner_lons", Float64, ("XCdim", "YCdim", "nf"),
                      attrib=["units" => "degrees_east", "long_name" => "longitude"])
    ds_clats = defVar(ds, "corner_lats", Float64, ("XCdim", "YCdim", "nf"),
                      attrib=["units" => "degrees_north", "long_name" => "latitude"])
    ds_clons[:, :, :] = Float64.(spec.corners_lon)
    ds_clats[:, :, :] = Float64.(spec.corners_lat)

    # CF grid-mapping anchor (data variables carry `grid_mapping = "cubed_sphere"`).
    defVar(ds, "cubed_sphere", Int32, ();
           attrib=["grid_mapping_name" => "gnomonic cubed-sphere"])

    # Global attributes: CF + GCHP-style plus self-describing CS provenance.
    ds.attrib["Conventions"] = "CF"
    ds.attrib["grid_mapping_name"] = "gnomonic cubed-sphere"
    ds.attrib["grid"] = "cubed-sphere"
    ds.attrib["Nc"] = spec.Nc
    ds.attrib["cs_coordinate_law"] = coordinate_law_tag(coordinate_law(spec.definition))
    ds.attrib["cs_center_law"] = center_law_tag(center_law(spec.definition))
    ds.attrib["cs_longitude_offset_deg"] = longitude_offset_deg(spec.definition)
    return nothing
end
