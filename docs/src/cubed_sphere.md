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

A grid is named `C<Nc>`, where `Nc` is the number of cells per panel edge; the
cube has `6·Nc²` cells total. For example, `Nc = 360` is **C360** (~25 km).

!!! note "Projection provenance"
    The closed-form gnomonic projection and its inverse (`lonlat_to_panel_xy`)
    are ported verbatim from the AtmosTransportModel package
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

Because the lon/lat → cube map is nonlinear and discontinuous across panel
edges, footprints are subdivided in **lon/lat space** (not in fractional-index
space as for rectangular grids). Each footprint corner quadrilateral is split
into `n × n` sub-points, each projected independently onto its `(panel, cell)`
via `lonlat_to_panel_xy`, contributing weight `1/n²`. A footprint straddling a
cube edge therefore correctly deposits weight onto both panels. Corner
longitudes are unwrapped to a common branch first so the subdivision is sane for
near-pole footprints that cross the ±180° seam.

The default oversampling factor on the cube is a small fixed value (`n = 4`);
override it with `--nOversample` / `SubpixelGridding(n_oversample=...)`.

!!! warning "Sequential CPU only"
    Cubed-sphere gridding currently runs on the sequential CPU path. The
    KernelAbstractions GPU kernels are 2D-grid-specific; a CS GPU path is a
    planned follow-up. Passing `--backend cuda/metal/cpu` with `--gridType cs`
    raises a clear error.

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
| Counts | `n(time, nf, Ydim, Xdim)` — pixel count per cell |

In xarray you can plot a face directly on its 2D coordinates:

```python
import xarray as xr
ds = xr.open_dataset("tropomi_sif_c360.nc")
ds["sif_743"].isel(time=0, nf=2).plot(x="lons", y="lats")   # one panel
```
