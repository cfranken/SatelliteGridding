# ---------------------------------------------------------------------------
# Output-grid geometry abstraction.
#
# `AbstractGridSpec` is the dispatch surface that lets `gridder.jl`,
# `ncio.jl`, and `oversampling.jl` serve multiple output grid geometries
# (regular lat/lon, cubed sphere, reduced Gaussian, …) without case
# analysis. PR #1 introduces the surface and a `RectangularGridSpec`
# implementation that mirrors the existing inlined logic; subsequent PRs
# add `CubedSphereGridSpec` and rewire the call sites through these
# methods.
#
# Kept separate from `types.jl` so the configuration / metadata records
# (DataSourceConfig, TimeSpec, FilterRule, AbstractGriddingMethod) stay
# independent of the grid-extensibility contract — mirrors the
# `AbstractMeshes.jl` + per-mesh-type file layout used by
# AtmosTransport.
# ---------------------------------------------------------------------------

"""
    AbstractGridSpec{T<:AbstractFloat}

Common supertype for output-grid specifications. Subtypes describe how a
satellite observation's `(lat, lon)` maps onto cell indices on a
particular grid geometry and how the matching output NetCDF file is laid
out.

Every concrete subtype `G <: AbstractGridSpec{T}` must implement five
methods so the same gridder/oversampling/ncio machinery serves all grids
via multiple dispatch:

| Method | Returns | Purpose |
|--------|---------|---------|
| [`to_fractional_index`](@ref)`(spec, lat, lon)` | `Tuple` of fractional indices | Project a point onto continuous cell indices |
| [`in_bounds`](@ref)`(spec, lat, lon)` | `Bool` | Cheap scalar bounding-box pre-filter |
| [`grid_shape`](@ref)`(spec)` | `NTuple{N, Int}` | Spatial size of the output accumulator array |
| [`payload_dims`](@ref)`(spec)` | `NTuple{N, String}` | NetCDF dimension names for payload variables |
| [`write_coordinates!`](@ref)`(ds, spec)` | `Nothing` | Write the spatial coordinate axes / aux coords to a NetCDF dataset |

`N` is the spatial dimensionality of the accumulator: `2` for a regular
lat/lon grid (`(n_lon, n_lat)`), `3` for a cubed sphere (`(Nc, Nc, 6)`).
The fractional-index tuple is allowed to mix types when an axis is
categorical (e.g. cubed sphere returns `(panel::Int, s::T, t::T)`),
so the return type is documented as `Tuple` rather than `NTuple{N, T}`.
The first `N` integer cell indices recovered via `floor(Int, …)` along
spatial axes match `grid_shape(spec)` axis-for-axis.
"""
abstract type AbstractGridSpec{T<:AbstractFloat} end

# Interface stubs. Concrete subtypes add methods; an unspecialized call
# throws a clear `MethodError` rather than falling back silently.

"""
    to_fractional_index(spec::AbstractGridSpec, lat, lon) -> Tuple

Map a single `(lat, lon)` in degrees to continuous fractional cell
indices on `spec`. The return tuple has the **same axis order as
[`grid_shape`](@ref)**, so the integer cell index along axis `k` is
`floor(Int, result[k])` and (when the call site has filtered with
[`in_bounds`](@ref)) lies in `1..grid_shape(spec)[k]`.

For [`RectangularGridSpec`](@ref) the return is `(ilon, ilat)` — **lon
first**, matching `grid_shape = (n_lon, n_lat)`. For a cubed-sphere
grid (PR #2) the return is `(s, t, panel::Int)` with `(Nc, Nc, 6)` as
the grid shape; an axis whose index is categorical is allowed to be
`Int` rather than the float type `T`.
"""
function to_fractional_index end

"""
    in_bounds(spec::AbstractGridSpec, lat, lon) -> Bool

Cheap scalar bounding-box pre-filter, intended to be broadcast at call
sites: `any(in_bounds.(spec, lat_centers, lon_centers))`. Returns
`true` when the point is strictly inside the grid's bounding box (open
intervals on both axes, matching the existing call-site convention).
Whole-sphere grids (cubed sphere, reduced Gaussian) return `true` for
every finite input.
"""
function in_bounds end

"""
    grid_shape(spec::AbstractGridSpec) -> NTuple{N, Int}

Spatial dimensions of the output accumulator array, in the same axis
order as [`to_fractional_index`](@ref).
"""
function grid_shape end

"""
    payload_dims(spec::AbstractGridSpec) -> NTuple{N, String}

NetCDF dimension names used when defining payload variables on this
grid. Returned in the same order as [`grid_shape`](@ref).
"""
function payload_dims end

"""
    accumulator_shape(spec::AbstractGridSpec) -> NTuple{M, Int}

Spatial dimensions of the in-memory accumulator array that the gridder
allocates, in the same axis order as the indices returned by
[`to_fractional_index`](@ref). This is the layout the Welford / scatter
machinery indexes with `(col, row)`.

For grids whose accumulator matches the NetCDF spatial layout this is just
[`grid_shape`](@ref) (the default). A cubed-sphere grid folds its
`(Nc, Nc, 6)` cube into a 2D `(Nc, 6·Nc)` accumulator and un-folds only at
write time, so its `accumulator_shape` equals its (already-folded)
`grid_shape`.
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
# These mirror the closed-form affine map currently inlined inside
# `_process_l2_file!` / `_process_l2_file_ka!` and the dim/coord-var
# definitions in `create_output_dataset`. The call sites in `gridder.jl`
# and `ncio.jl` will be rewired through these methods in a follow-up PR;
# this PR introduces the dispatch surface only, so existing performance
# and behavior are unchanged.

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
