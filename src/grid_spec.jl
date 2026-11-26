# ---------------------------------------------------------------------------
# Output-grid geometry abstraction.
#
# `AbstractGridSpec` is the dispatch surface that lets `gridder.jl`,
# `ncio.jl`, and `oversampling.jl` serve multiple output grid geometries
# (regular lat/lon, cubed sphere, HEALPix, reduced Gaussian, …) without
# case analysis. The rectangular implementation mirrors the original
# inlined affine map. Spherical grids locate individual subpoints.
#
# Kept separate from `types.jl` so the configuration / metadata records
# (DataSourceConfig, TimeSpec, FilterRule, AbstractGriddingMethod) stay
# independent of the grid-extensibility contract. This mirrors the
# `AbstractMeshes.jl` + per-mesh-type file layout used by AtmosTransport.
# ---------------------------------------------------------------------------

"""
    AbstractGridSpec{T<:AbstractFloat}

Common supertype for output-grid specifications. Subtypes describe how a
satellite observation's `(lat, lon)` maps onto cell indices on a
particular grid geometry and how the matching output NetCDF file is laid
out.

Concrete subtypes `G <: AbstractGridSpec{T}` supply their geometry to
the same gridder, oversampling, and NetCDF machinery via multiple dispatch.

| Method | Returns | Purpose |
|--------|---------|---------|
| [`to_fractional_index`](@ref)`(spec, lat, lon)` | `Tuple` of fractional indices | Project a point onto continuous cell indices where the geometry supports them |
| [`in_bounds`](@ref)`(spec, lat, lon)` | `Bool` | Cheap scalar bounding-box pre-filter |
| [`grid_shape`](@ref)`(spec)` | `NTuple{N, Int}` | Spatial grid size in index order, used by default for the accumulator |
| [`payload_dims`](@ref)`(spec)` | `NTuple{N, String}` | NetCDF dimension names for payload variables |
| [`write_coordinates!`](@ref)`(ds, spec)` | `Nothing` | Write the spatial coordinate axes / aux coords to a NetCDF dataset |

For grids with fractional indices, flooring the components gives one-based
integer indices in `grid_shape(spec)` axis order. The return type is
`Tuple` because an axis with categorical indices may use a different type.
The spatial dimensionality `N` of the accumulator is 2 for regular lat/lon
grids (`(n_lon, n_lat)`) and for cubed spheres (which fold `(Nc, Nc, 6)`
into `(Nc, 6·Nc)`).

Whole-sphere grids locate subpoints directly through the
[`AbstractSphericalGridSpec`](@ref) interface. Their on-disk dimensions
can differ from the folded accumulator shape.
"""
abstract type AbstractGridSpec{T<:AbstractFloat} end

"""
    AbstractSphericalGridSpec{T} <: AbstractGridSpec{T}

Whole-sphere grid for quadrilateral footprint gridding.
Implement these methods to use a new spherical grid with [`grid`](@ref).

- `_spherical_cell(spec, lat, lon)` locates a point in degrees and returns
  one-based `(col, row)` indices in the accumulator.
- `grid_shape(spec)` returns the two-dimensional accumulator shape, or
  implement `accumulator_shape(spec)` directly to specify a folded layout.
- `_output_spatial_shape(spec)` returns the NetCDF spatial shape in Julia
  dimension order. The writer reshapes the accumulator in linear index order.
- `payload_dims(spec)` returns the corresponding NetCDF spatial dimension names.
- `write_coordinates!(ds, spec)` defines those dimensions and writes coordinates.

Optionally specialize `_spherical_subdivision(spec)`. The default
`Val(:lonlat)` uses the corner shortcut and samples in longitude/latitude
space. `Val(:unit_vectors)` samples every footprint in 3-D and can be used
with a fixed `n_oversample` (four per side by default). HEALPix also provides
an adaptive sample count.

To opt into KernelAbstractions, define `_spherical_subdivision(spec) = Val(:unit_vectors)`
so the sequential and device paths both use unit-vector subdivision. Implement
`_device_locator(spec)` to return an isbits locator and
`locate(locator, x, y, z)` to return a one-based linear accumulator index
for a normalized unit vector. The locator must accept the
chosen geometry float type and use allocation-free device operations.
Optionally implement `_device_cell_spacing(locator, F)` to supply the angular
grid spacing for adaptive subdivision. The default uses four samples per side.
"""
abstract type AbstractSphericalGridSpec{T<:AbstractFloat} <: AbstractGridSpec{T} end

function _spherical_cell end
function locate end

_device_locator(::AbstractSphericalGridSpec) = nothing
_device_cell_spacing(locator, ::Type{F}) where {F} = F(Inf)

# All device radicands are constructed from finite, nonnegative terms.
@inline _spherical_sqrt(x::Union{Float32,Float64}) = Core.Intrinsics.sqrt_llvm(x)

_spherical_subdivision(::AbstractSphericalGridSpec) = Val(:lonlat)
_spherical_subpixel_type(::AbstractSphericalGridSpec) = Float64

@inline function _canonical_longitude(lon::Float64)
    wrapped = mod(lon, 360.0)
    # A negative angle smaller than half an ulp at 360° can round up to 360°.
    return wrapped == 360.0 ? 0.0 : wrapped
end

# Interface stubs. Concrete subtypes add methods; an unspecialized call
# throws a clear `MethodError` rather than falling back silently.

"""
    to_fractional_index(spec::AbstractGridSpec, lat, lon) -> Tuple

Map a single `(lat, lon)` in degrees to continuous fractional cell
indices on `spec`. The return tuple has the **same axis order as
[`grid_shape`](@ref)**, so the integer cell index along axis `k` is
`floor(Int, result[k])` and (when the call site has filtered with
[`in_bounds`](@ref)) lies in `1..grid_shape(spec)[k]`.

For [`RectangularGridSpec`](@ref) the return is `(ilon, ilat)`, with
**longitude first**, matching `grid_shape = (n_lon, n_lat)`. For a
cubed-sphere grid the return is `(col, row)` in the folded `(Nc, 6·Nc)`
accumulator. A categorical index may be an `Int` alongside fractional
indices of type `T`. Grids with direct cell locators such as HEALPix use
`_spherical_cell` for footprint gridding.
"""
function to_fractional_index end

"""
    in_bounds(spec::AbstractGridSpec, lat, lon) -> Bool

Cheap scalar bounding-box pre-filter, intended to be broadcast at call
sites: `any(in_bounds.(spec, lat_centers, lon_centers))`. Returns
`true` when the point is strictly inside the grid's bounding box (open
intervals on both axes, matching the existing call-site convention).
Whole-sphere grids have periodic longitude. The cubed-sphere implementation
checks finiteness, and HEALPix also checks the latitude range here. Spherical
footprint gridding validates the centre and all four vertices together
during input preparation.
"""
function in_bounds end

"""
    grid_shape(spec::AbstractGridSpec) -> NTuple{N, Int}

Spatial grid dimensions in index order. For grids with
[`to_fractional_index`](@ref), the axis order matches the fractional
indices. The gridder allocates this shape by default, and a grid may
specialize [`accumulator_shape`](@ref) separately.
"""
function grid_shape end

"""
    payload_dims(spec::AbstractGridSpec) -> NTuple{N, String}

NetCDF dimension names used when defining payload variables on this grid.
The spatial dimensions are returned in Julia dimension order, matching
`_output_spatial_shape(spec)`. Rectangular grids use the same order as
[`grid_shape`](@ref), and the writer appends time after the spatial axes.
"""
function payload_dims end

"""
    accumulator_shape(spec::AbstractGridSpec) -> NTuple{M, Int}

Spatial dimensions of the in-memory accumulator array that the gridder
allocates, in the same axis order as the indices returned by
[`to_fractional_index`](@ref) or `_spherical_cell`. The Welford / scatter
machinery indexes this layout with `(col, row)`.

The default is [`grid_shape`](@ref). A cubed-sphere grid folds its
`(Nc, Nc, 6)` cube into a 2D `(Nc, 6·Nc)` accumulator and unfolds only
at write time. HEALPix folds the NESTED sequence into `(nside^2, 12)`
in memory and writes it on a single cell axis. Both use the default
because their `grid_shape` is already folded.
"""
accumulator_shape(spec::AbstractGridSpec) = grid_shape(spec)

"""
    write_coordinates!(ds, spec::AbstractGridSpec) -> Nothing

Define the spatial dimensions and write the coordinate axis variables
(plus auxiliary `lons`/`lats` arrays, `corner_lons`/`corner_lats`, and a
`grid_mapping` scalar for grids that need them) onto an open NetCDF
dataset `ds`. Time and payload variables are defined separately by the
caller.
"""
function write_coordinates! end

"""
    RectangularGridSpec{T<:AbstractFloat} <: AbstractGridSpec{T}

Output grid on a regular latitude/longitude rectangle with uniform cell
spacing. The legacy name `GridSpec` is a `const` alias of this type for
back-compatibility.

# Fields
- `lat_min`, `lat_max`: Latitude bounds (degrees)
- `lon_min`, `lon_max`: Longitude bounds (degrees)
- `dlat`, `dlon`: Grid cell size (degrees)
- `lat`, `lon`: Vectors of cell center coordinates
"""
struct RectangularGridSpec{T<:AbstractFloat} <: AbstractGridSpec{T}
    lat_min::T
    lat_max::T
    lon_min::T
    lon_max::T
    dlat::T
    dlon::T
    lat::Vector{T}
    lon::Vector{T}
end

"""
    RectangularGridSpec(; lat_min=-90f0, lat_max=90f0, lon_min=-180f0,
                          lon_max=180f0, dlat=1f0, dlon=1f0)

Construct a `RectangularGridSpec` with cell centers computed from bounds
and resolution. `GridSpec(...)` is accepted as a synonym.
"""
function RectangularGridSpec(; lat_min::T=-90.0f0, lat_max::T=90.0f0,
                              lon_min::T=-180.0f0, lon_max::T=180.0f0,
                              dlat::T=1.0f0, dlon::T=1.0f0) where {T<:AbstractFloat}
    eps = dlat / 100
    lat = collect(lat_min + dlat / 2:dlat:lat_max - dlat / 2 + eps)
    lon = collect(lon_min + dlon / 2:dlon:lon_max - dlon / 2 + eps)
    RectangularGridSpec{T}(lat_min, lat_max, lon_min, lon_max, dlat, dlon, lat, lon)
end

"Back-compat alias for [`RectangularGridSpec`](@ref). All call sites using `GridSpec(...)` and `GridSpec{T}` keep working unchanged."
const GridSpec = RectangularGridSpec

# ---------------------------------------------------------------------------
# AbstractGridSpec interface — RectangularGridSpec implementations
# ---------------------------------------------------------------------------
#
# These mirror the original closed-form affine map in
# `_process_l2_file!` / `_process_l2_file_ka!` and the dimension and
# coordinate-variable definitions in `create_output_dataset`. The
# rectangular L2 paths use the affine map in bulk, while scalar callers
# use this interface. The NetCDF writer dispatches through the methods below.

@inline function to_fractional_index(spec::RectangularGridSpec{T},
                                     lat::Real, lon::Real) where {T}
    n_lat = length(spec.lat)
    n_lon = length(spec.lon)
    ilon = T((lon - spec.lon_min) / (spec.lon_max - spec.lon_min) * n_lon) + one(T)
    ilat = T((lat - spec.lat_min) / (spec.lat_max - spec.lat_min) * n_lat) + one(T)
    return (ilon, ilat)
end

@inline function in_bounds(spec::RectangularGridSpec, lat::Real, lon::Real)
    return spec.lat_min < lat < spec.lat_max &&
           spec.lon_min < lon < spec.lon_max
end

grid_shape(spec::RectangularGridSpec) = (length(spec.lon), length(spec.lat))

payload_dims(::RectangularGridSpec) = ("lon", "lat")

function write_coordinates!(ds, spec::RectangularGridSpec)
    defDim(ds, "lon", length(spec.lon))
    defDim(ds, "lat", length(spec.lat))
    ds_lat = defVar(ds, "lat", Float32, ("lat",),
                    attrib=["units" => "degrees_north", "long_name" => "Latitude"])
    ds_lon = defVar(ds, "lon", Float32, ("lon",),
                    attrib=["units" => "degrees_east", "long_name" => "Longitude"])
    ds_lat[:] = spec.lat
    ds_lon[:] = spec.lon
    return nothing
end
