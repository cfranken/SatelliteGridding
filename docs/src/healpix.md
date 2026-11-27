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
per-footprint cell overlaps also supply the optional `sum_w2` and
`n_soundings` statistics.

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
shifted sums, deviations, and cached daily moments. Accumulators default to
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
shifted first pass. Float32 KA runs replay the files to refine each cell's
mean using deviations around it. Float64 KA runs replay when standard
deviation is requested. [GPU Acceleration](@ref) describes precision limits,
memory use, and the spherical grid locator interface.

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
standard deviation when the moments are finite. Spherical
KA writes fill when the variance or its inputs become non-finite in
accumulator precision. The calculation uses no Bessel correction.
With Float32 accumulators, including on Metal, roundoff can erase small
standard deviations, especially in cells with many contributions.
See [the GPU precision discussion](gpu.md#accumulation-and-standard-deviation).

NetCDF stores `n = 0` for unobserved cells, including entirely empty time
slices, and reserves `-1` for unwritten weights. A cell with `n > 0` can
still contain science fill, for example when its standard deviation overflows.

L2 NetCDF records `accumulator_precision` and `geometry_precision` as global
attributes. For gridding inputs to coarsening or regridding, use the
[exact-archive recipe](#exact-archive-recipe).

With the optional `sum_w2` and `n_soundings` fields, fractional overlaps act
as reliability weights for independent soundings with a common mean and variance.
Require `sum_w2 > 0`, `n_soundings ≥ 2`, and `N_eff > 1`, where Kish
`N_eff = n² / sum_w2`. The unbiased sample variance is
`s² = M2 / (n − sum_w2/n) = std² * N_eff / (N_eff − 1)` and the standard
error is `SE = std / √(N_eff − 1)`. Avoid `M2 / (n − 1)` and `std / √n`.
Overlapping rolling slices reuse soundings, and soundings shared across
cells induce covariance. [Overlap statistics configuration](config.md#netcdf-storage-and-provenance).

All rolling windows with standard deviation merge centered daily moments
in Float64, as do mean-only Float32 KA runs. Float64 KA restores the shared
reference after merging. Float32 KA replay supplies the bounded absolute
mean to the merge.

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

## Coarsening a NESTED archive

Use `coarsen_healpix("fine.nc", "coarse.nc"; levels=2)` in Julia, or
pass `nside=…` to select a parent resolution directly. Each level halves
`nside` and combines four consecutive NESTED cells into one parent.
The output includes every time slice, the original time and window
coordinates, and recomputed HEALPix coordinates and bounds (if present).
The global attributes `coarsen_source`, `coarsen_source_nside`,
`coarsen_target_nside`, and `coarsen_levels` record the original archive,
its resolution, the output resolution, and the total levels removed across
all coarsening steps.

The coarsener uses `satellitegridding_role` to identify weights and science
fields and `satellitegridding_mean` and `satellitegridding_std` to link each
mean and its deviation. These tags take precedence over `legacy_std`.
Untagged archives without `<name>` and `<name>_std` pairs coarsen as means.
When a pair exists, the default `legacy_std=:error` stops before creating
output because `<name>_std` could be an independent science field.

- Use `legacy_std=:by_name` for archives written with `compute_std=true`.
  Every science mean must have a matching population standard deviation
  on `(cell, time)`, and the pairs must be disjoint.
- Use `legacy_std=:none` to coarsen every science field as an independent
  mean, including fields whose names end in `_std`.

For each parent, the coarsener sums child weights and combines child means
and central moments with Chan's formula in Float64. It recovers each
child's central moment as `M2 = n * std^2` and writes the parent's weighted
population standard deviation as `sqrt(M2 / n)`. It skips cells with
`n = 0`. A missing mean in a positive-weight child propagates to the parent
mean and standard deviation for that variable. A missing child standard
deviation propagates to the parent standard deviation.

### Exact-archive recipe

Grid the fine archive with Float64 accumulation and spherical geometry.
In Julia, use `HEALPixGridSpec(; nside=..., T=Float64)` for sequential
gridding. On the sequential CLI, use `--accumulator-precision float64`.
KA CPU and CUDA default to Float64 for both precisions. Enable standard
deviations with `compute_std=true` in Julia or `--compSTD` on the CLI.

```toml
[output]
precision = "float64"
significant_digits_mean = "none"
significant_digits_std = "none"
fill_value = nan
```

Coarsening preserves the moments of the archived sub-sample set up to
Float64 summation and serialization roundoff, including with adaptive
oversampling. To reproduce those moments by gridding directly at the parent
resolution, use the same fixed `n_oversample` (`--nOversample` on the CLI),
geometry precision, and time windows for the two runs. Adaptive oversampling
can choose different sub-samples at each resolution. Quantization of nearly
equal serialized child means can dominate a small parent standard deviation.
Coarsening writes Float64 without decimal rounding and warns about Float32
accumulation, Float32 storage, decimal rounding, and missing precision or
rounding metadata.

### Coarsened overlap statistics

When present, `sum_w2` and `n_soundings` are summed across descendants and
marked with `coarsening_approximation`. A missing statistic in any observed
child propagates to that parent statistic. For `L` total levels removed,
compare with parent statistics for the same archived sub-samples. The summed
`sum_w2` is between `4⁻ᴸ` and `1` times the true parent value, while the summed
`n_soundings` is between `1` and `4ᴸ` times the true parent count. A footprint
touching multiple descendants contributes multiple counts and loses the
cross-products of its weights. This can inflate `N_eff` by up to `4ᴸ`.
**Do not use these sums for sample variance or SE.** Even one footprint can
pass both estimator guards after coarsening and give a spurious zero SE.
Grid directly at the parent resolution for these estimates.

In a one-day OCO-2 run on CUDA with Float64 and no decimal rounding, using
`n_oversample=4`, coarsening from `nside=512` to `128` matched direct
`nside=128` gridding with zero difference in `n`, a maximum mean difference
of `6 × 10⁻¹⁴ ppm`, and a maximum relative standard deviation difference of
`9 × 10⁻¹⁴`. The summed `n_soundings` was 20% higher (97,157 vs 81,077) and
the summed `sum_w2` was 6.6% lower (69,099 vs 73,950).

## Known limitations and planned work

HEALPix supports quadrilateral Level-2 footprints to NetCDF. Circular
footprints, centre-coordinate gridding such as MODIS, Zarr output, and
UGRID export are outside this support.

- Per-cell WGS84 ellipsoidal areas are planned separately.
