# HEALPix Grids

[`HEALPixGridSpec`](@ref) provides global gridding in HEALPix NESTED order
with equal solid angle per cell on the sphere. The resolution is
`nside = 2^level`, with `12 * nside^2` cells and zero-based Int64 cell IDs.
The four children of cell `p` at the next level have IDs
`4p:4p+3`. The parent ID of a child `c` is `c >> 2`.
[HEALPix geometry and numbering, Górski et al. (2005)](https://arxiv.org/abs/astro-ph/0409513).

Input geodetic latitudes are used directly as spherical latitudes, as for the
cubed-sphere grid, so the grids use the same latitude convention. True WGS84
cell areas vary by about 1.35% from equator to pole, from −0.45% at the
equator to +0.90% at the poles relative to the spherical `cell_area`.
Longitude is periodic. Input longitude is canonicalized in degrees before
locating cells, and cell centre and vertex longitudes are in `[0, 360)`.

## Usage

```julia
using SatelliteGridding, Dates

config = load_config("examples/tropomi_sif.toml")
spec = HEALPixGridSpec(nside=64, cell_bounds=true)
times = TimeSpec(DateTime("2020-07-01"), DateTime("2020-07-16"), Day(16))
grid_l2(config, spec, times; n_oversample=8, compute_std=true,
        outfile="tropomi_healpix.nc")
```

```bash
julia --project=. bin/grid.jl l2 \
    --config examples/tropomi_sif.toml \
    --gridType healpix --nside 64 --nOversample 8 \
    --startDate 2020-07-01 --stopDate 2020-07-16 --dDays 16 \
    --outFile tropomi_healpix.nc
```

Every footprint is sampled, including those whose four corners fall in a
single cell. Corners are converted to unit vectors and ordered counterclockwise
in the tangent plane at their normalized vector sum. The four vectors are
interpolated bilinearly at `((i - 0.5)/n, (j - 0.5)/n)` and the resulting
vectors are normalized to locate the subpoints on the sphere. This handles
footprints crossing the dateline and those enclosing a pole.

Subpoints carry equal weight `1/n²` in the bilinear unit-vector
parameterization, corresponding to uniform angular sampling within the
detector's instantaneous field of view (IFOV). For parallelogram footprints
this gives spherical-area weighting to footprint-scale accuracy, with a
small correction from radial normalization. For tapered footprints the
difference follows the taper. For realistic OCO and TROPOMI pushbroom taper
of up to about 0.5%, the estimated difference in cell weights is about 0.1%
or less, roughly 30 times smaller than finite-sampling error in the swath-edge
tests.

The gridding loop counts the `k` subpoints in each cell and updates its
statistics once per footprint with weight `k/n²`. Footprints entirely
inside one cell contribute exactly one unit of weight. The resulting
per-footprint cell overlaps also provide the weights needed for a future
sum of squared weights (`Σw²`).

By default `n = clamp(ceil(Int, 4 * L_max / d_cell), 4, 24)` per footprint.
`L_max` is the largest great-circle distance between any two footprint
corners and `d_cell = sqrt(4πR² / npix)` is the mean cell spacing on a sphere
of radius `R`. An explicit `n_oversample` in Julia or `--nOversample` on the
CLI uses a fixed integer from 1 to 1024 samples per side.
The upper clamp can leave interior cells unsampled. Use an explicit
`n_oversample` for coarse sounders with `L_max ≳ 75 km` at `nside = 2048`,
or for TROPOMI swath-edge footprints at `nside ≥ 8192`.

## GPU and KernelAbstractions

Use `backend=resolve_backend("cuda")`, `resolve_backend("metal")`, or
`resolve_backend("cpu")` in Julia, or `--backend cuda`, `--backend metal`,
or `--backend cpu` on the CLI. All three support quadrilateral footprints,
tiled and rolling time windows, and `compute_std=true` (`--compSTD`).
Standard deviation rereads the files in each tile or rolling-cache day.

CUDA and KA CPU use `Float64` geometry by default to match sequential cell
assignment. Set `geometry_precision=Float32` (`--geometry-precision float32`)
to use single precision. Set `accumulator_precision`
(`--accumulator-precision float32` or `float64`) independently for weights,
shifted sums, deviations, and rolling moments. Accumulators default to
Float64 on CUDA and KA CPU. Metal uses Float32 for both settings and
rejects Float64.
Float32 geometry quantizes the input corner vectors and orders them in the tangent
plane in the same precision. It can move samples across cell edges or
change the adaptive sample count near its rounding thresholds.
At `nside=2048`, the Float32 cell-edge shift is at most about 0.33 m.
For representative OCO and TROPOMI footprints, the total subpoint
displacement is about 0.2 m RMS and at most about 1.1 m, below the typical
50–500 m L2 geolocation uncertainty. In sparsely grazed perimeter cells (`n ≪ 1`), one flipped
subpoint can change occupancy or shift the mean by an amount comparable
to the scene standard deviation. In the XCO₂ examples, cell means with
`n ≳ 2` shifted by about 0.05 ppm or less.
The first nonempty batch supplies a per-variable median reference for the
shifted first pass. Mean-only runs skip refinement. In either accumulator
precision `A`, values with magnitude below about
`eps(A) * abs(reference - value)` can round away in that first pass.
`eps(Float64)` is about `2e-16`, and `eps(Float32)` is about `1e-7`.
Enable `compute_std` to recover those values by refining each cell's mean
using deviations around it. [GPU Acceleration](@ref) gives examples of
precision limits and describes memory use and the spherical grid locator interface.

## NetCDF layout

Gridded variables, weight sums `n`, and optional `*_std` all use
`(cell, time)` in Julia and `(time, cell)` in xarray or CDL. They reference the
`healpix` grid mapping and the `lat lon cell_ids` coordinates. Files declare
`Conventions = "CF-1.13"`. The `cell_ids` variable holds the NESTED IDs with
`standard_name = "healpix_index"`, `units = "1"`, and xdggs metadata
(`grid_name = "healpix"`, `level`, `indexing_scheme = "nested"`, and
`ellipsoid = "sphere"`). The scalar `cell_area` gives the common area in square
metres on a sphere of radius 6,371 km.
[CF 1.13, Appendix F](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.13/cf-conventions.html#_healpix).

`n` is dimensionless and sums fractional footprint overlap weights, so a
sounding split across cells contributes a fraction to each cell. The optional
`*_std` is the weighted population standard deviation, `sqrt(M2 / n)`, where
`M2` is the weighted sum of squared deviations from the cell mean. It has the
source variable's units. A cell touched by a single footprint has zero
standard deviation in tiled mode when the moments are finite. Spherical
KA writes fill when the variance or its inputs become non-finite in
accumulator precision. The calculation uses no Bessel correction.
With Float32 accumulators, including on Metal, a standard deviation below
roughly one ulp of the cell mean (`σ/|mean| ≲ 1e-7`) can be reported as zero.
Use Float64 accumulators (the CUDA and KA CPU default) when such small
relative variability matters.

NetCDF stores `n = 0` for unobserved cells, including entirely empty time
slices, and reserves `-1` for unwritten weights. A cell with `n > 0` can
still contain science fill, for example when its standard deviation overflows.

To archive Float64 results exactly, accumulate in Float64 and set
`[output] precision = "float64"`, `significant_digits_mean = "none"`, and
`significant_digits_std = "none"`. In Julia, use
`HEALPixGridSpec(; nside=..., T=Float64)` for sequential gridding. The CLI
creates Float32 grids.
HEALPix on `--backend cpu` or `--backend cuda` accumulates in Float64 by
default. Other CLI modes accumulate in Float32. Use this Float64 configuration
for outputs that feed HEALPix 4:1 coarsening or regridding. At fine `nside`,
quantization of nearly equal child means can dominate a small parent standard
deviation.

Spherical KA rolling windows combine centered daily moments using corrected
means relative to the run reference. Sequential
rolling windows reconstruct variance from unshifted raw moments, which
can lose small variances in Float32 when the absolute tracer values are
large. Use `HEALPixGridSpec(nside=...; T=Float64)` for sequential rolling
standard deviations of high-offset tracers such as XCO₂.

For xdggs, open the file with xarray and use
`xdggs.decode(ds.set_coords("cell_ids"))` to decode the ID coordinate.
[xdggs HEALPix tutorial](https://xdggs.readthedocs.io/en/latest/tutorials/healpix.html).

Float64 `lon` and `lat` hold cell centres. Set `cell_bounds=true` in Julia to
also write Float64 `lon_bounds` and `lat_bounds` with Julia dimensions
`(nv, cell)`. The four corner vertices follow healpy's north, west, south,
east order and are intended for plotting. Conservative regridding needs the
curved HEALPix edges. Connecting the vertices by great circles changes the
area of cells touching a pole by about 10%. Coordinate arrays are generated
in chunks during file creation.
Gridding uses dense accumulators, so memory still grows with the cell count
and the number of variables.

## Known limitations and planned work

HEALPix supports quadrilateral Level-2 footprints to NetCDF. Circular
footprints, centre-coordinate gridding such as MODIS, Zarr output, UGRID
export, and coarsening or aggregation helpers are outside this support.

- In sequential rolling windows, Float32 raw moment cancellation can
  erase variance or introduce spurious variance for high-offset tracers.
- Outputs contain neither per-cell `sum_w2` nor counts of contributing soundings.
- Per-cell WGS84 ellipsoidal areas are planned separately.
- Cubed-sphere seam corner ordering, pole-enclosing footprints, and
  longitude/latitude subdivision are planned for a separate review.
