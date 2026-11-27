# Cubed-Sphere Grids (GEOS / GCHP)

In addition to regular latitude/longitude grids ([`RectangularGridSpec`](@ref)),
SatelliteGridding can grid Level-2 footprints directly onto a **cubed-sphere**
output grid in the GEOS / GCHP convention via [`CubedSphereGridSpec`](@ref). This
is useful for comparing satellite observations against GEOS-Chem (GCHP) output on
its native grid without an intermediate regridding step.

## Convention

The default geometry is [`GMAOCubedSphereDefinition`](@ref) — the grid used by
GEOS-IT / GEOS-FP and by GCHP when driven by GMAO meteorology:

- **equal-distance gnomonic** edge coordinate law (`grid_type = 0` in GEOS/FV),
- **four-corner-normalized** cell centers (GEOS `cell_center2`),
- **GEOS-native panel ordering** (panels 1–2 equatorial, 3 north pole, 4–5
  equatorial, 6 south pole), and
- a final **−10° longitude shift**.

An [`EquiangularCubedSphereDefinition`](@ref) (classic tan-spaced equiangular
gnomonic, no offset) is also available for synthetic targets.

A grid is named `C<Nc>`, where `Nc` is the number of cells per panel edge.
The cube has `6·Nc²` cells total. For example, `Nc = 360` is **C360** (~25 km).
Input geodetic latitudes are used directly as spherical latitudes, as on
the HEALPix grid.

!!! note "Projection provenance"
    The closed-form gnomonic projection and its inverse (`lonlat_to_panel_xy`)
    are adapted from the AtmosTransportModel package
    (`src/Grids/CubedSphereMesh.jl`, commit `6522d85`). See
    [`src/cubed_sphere.jl`](https://github.com/cfranken/SatelliteGridding/blob/main/src/cubed_sphere.jl)
    for the full attribution.

## Library usage

```julia
using SatelliteGridding, Dates

config    = load_config("examples/tropomi_sif.toml")
grid_spec = CubedSphereGridSpec(; Nc = 360)             # C360, GMAO/GCHP default
time_spec = TimeSpec(DateTime("2020-07-01"), DateTime("2020-07-16"), Day(16))

grid(config, grid_spec, time_spec, SubpixelGridding();
     outfile = "tropomi_sif_c360.nc", keep_going = true)
```

The data-source config (`tropomi_sif.toml`) is grid-agnostic — only the
`grid_spec` changes between a regular and a cubed-sphere run. See
`examples/tropomi_sif_cubedsphere.jl` for a runnable C24 demo.

## CLI usage

```bash
julia --project=. bin/grid.jl l2 \
    --config examples/tropomi_sif.toml \
    --gridType cs --Nc 360 \
    --startDate 2020-07-01 --stopDate 2020-07-16 --dDays 16 \
    -o tropomi_sif_c360.nc
```

`--csConvention` selects `gmao` (default) or `equiangular`. For long global runs,
`bin/run_tropomi_c360.sh` fans the work out across time chunks (one Julia process
per chunk) and concatenates with `ncrcat`:

```bash
NC=360 JOBS=8 bash bin/run_tropomi_c360.sh
```

## How footprints are gridded

The cube and [HEALPix](healpix.md) share the same spherical footprint
subdivision. The four corners become unit vectors and are ordered in the
footprint's tangent plane. The sampler evaluates the centres of an `n × n`
bilinear parameter grid and normalizes each interpolated vector. The cube
locator projects each sample independently, so the same subdivision works
across the antimeridian and panel edges and for footprints that enclose a pole.
A cell hit by `k` samples receives one update with weight `k/n²`, its sampled
fraction of the footprint's bilinear parameter domain. The weights for each
valid footprint sum to one. Equal subpoint weights correspond to uniform angular
sampling within the detector's instantaneous field of view (IFOV). For realistic
OCO and TROPOMI pushbroom taper of up to about 0.5%, the estimated difference from
spherical-area cell weights is about 0.1% or less. See the [HEALPix weighting
discussion](healpix.md) for details.

Automatic subdivision chooses `n = clamp(ceil(4 L_max / Δ), 4, 24)` for each
footprint. Here `L_max` is the largest great-circle distance between any pair of
corners and `Δ = sqrt(4πR²/(6 Nc²))` is the mean cell spacing on a sphere of
radius `R`. At high resolution the largest cells have about `4/√3 ≈ 2.31`
times the area of the smallest on GMAO and `√2 ≈ 1.41` times on equiangular
grids. The narrowest GMAO corner cells are about `0.74 Δ` wide. Below the cap,
the nominal sample spacing is at most `Δ/4` along the longest span, giving
about three samples per narrow corner-cell width. The cap limits `n` to 24
once `L_max` exceeds `6 Δ`, about 154 km at C360 or 77 km at C720 on a sphere of
radius 6,371 km. Use `--nOversample` or `SubpixelGridding(n_oversample=...)` to
set `n` explicitly for wide footprints.

## Backends and precision

Quadrilateral gridding supports the sequential backend and
[KernelAbstractions](gpu.md) on KA CPU, CUDA, and Metal, including tiled and
rolling windows and standard deviation. Select an accelerator with
`backend=resolve_backend("cuda")` or `--backend cuda`.

Sequential geometry uses Float64. KA CPU and CUDA default to Float64 for both
geometry and accumulation. Set `geometry_precision` and `accumulator_precision`
independently in `grid_l2` (CLI `--geometry-precision` and
`--accumulator-precision`). Metal requires Float32 for both. As on HEALPix,
Float32 geometry rounds the corner vectors, interpolation and projection.
In a C720 test with 500,000 uniformly distributed points per definition, the
Float32 locator changed 0.0030% of GMAO cell assignments and 0.0034% of
equiangular assignments relative to Float64. A C720 scan of every panel's
horizontal and vertical grid lines at each cell-edge midpoint found assignment
changes up to about 1.42 m from a GMAO edge and 0.93 m from an equiangular edge,
including Float32 rounding of the input unit vectors. The GMAO inverse adds
opposite-signed terms near one panel edge to compute its fractional index, which
contributes to the shift. Exact edge and corner ties can choose a different
incident cell. Footprint interpolation and adaptive sample-count rounding can
also affect weights close to cell boundaries. At very large `Nc` (around `2^23`
and above), Float32 fractional cell indices can round by half a cell or more.
Use Float64 geometry at those resolutions.

The grid's `T` controls the precision of the precomputed cell centre and corner
coordinates. `CubedSphereGridSpec` with `T=Float32` computes these in Float32
and widens them to Float64 for NetCDF. Set `T=Float64` for double-precision
coordinates and for Float64 accumulation in tiled sequential runs.

## Output layout

The NetCDF output matches native GEOS / GCHP files, so it opens directly in
Panoply, xarray, and GCHP tooling:

| Item | Details |
|------|---------|
| Dimensions | `Xdim = Ydim = Nc`, `nf = 6`, `XCdim = YCdim = Nc+1`, unlimited `time` |
| Data variables | `(time, nf, Ydim, Xdim)`, with `coordinates = "lons lats"` and `grid_mapping = "cubed_sphere"` |
| Cell centers | 2D `lons(nf, Ydim, Xdim)`, `lats(...)` in degrees |
| Cell corners | `corner_lons(nf, YCdim, XCdim)`, `corner_lats(...)` |
| Fake 1D axes | `Xdim`/`Ydim` (degrees, GrADS compatibility); `nf = 1..6` |
| Weights | `n(time, nf, Ydim, Xdim)` stores the sum of fractional footprint weights |

In xarray you can plot a face directly on its 2D coordinates:

```python
import xarray as xr
ds = xr.open_dataset("tropomi_sif_c360.nc")
ds["sif_743"].isel(time=0, nf=2).plot(x="lons", y="lats")   # one panel
```
