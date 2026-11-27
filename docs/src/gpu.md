# GPU Acceleration

`SatelliteGridding.jl` supports accelerator execution through
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl).
The same gridding kernels can run on the KA CPU backend, NVIDIA GPUs through
[CUDA.jl](https://github.com/JuliaGPU/CUDA.jl), and Apple GPUs through
[Metal.jl](https://github.com/JuliaGPU/Metal.jl).

## Setup

CUDA.jl and Metal.jl are weak dependencies. They are loaded only when you
request that backend or explicitly import the package:

```julia
using SatelliteGridding

cuda_backend = resolve_backend("cuda")    # requires CUDA.jl
metal_backend = resolve_backend("metal")  # requires Metal.jl on macOS/Apple GPU
```

You can also use vendor packages directly:

```julia
using CUDA
using Metal

CUDA.devices()
Metal.devices()
```

## Usage

### Julia API

```julia
grid_l2(config, grid_spec, time_spec;
        backend = resolve_backend("cuda"),
        outfile = "output.nc")

grid_l2(config, grid_spec, time_spec;
        backend = resolve_backend("metal"),
        outfile = "output_metal.nc")
```

### CLI

```bash
julia --project=. bin/grid.jl l2 \
    --config examples/tropomi_sif.toml \
    --backend cuda \
    -o output.nc

julia --project=. bin/grid.jl l2 \
    --config examples/tropomi_sif.toml \
    --backend metal \
    -o output_metal.nc
```

Rectangular grids support quadrilateral footprints (`--footprint quad`) and circular
footprints (`--footprint circle`, `CircularFootprintGridding`). The circular path
uses a separate index kernel that masks samples outside the inferred circle or
ellipse before scatter accumulation. Both footprints support `compute_std=true`
(`--compSTD`) on KA CPU, CUDA, and Metal in tiled and rolling runs.

For rectangular, cubed-sphere, and HEALPix grids, set `accumulator_precision=Float32`
or `--accumulator-precision float32` to select single precision for weights,
shifted sums, and daily moments. KA CPU and CUDA default to Float64. Metal
requires Float32 and rejects an explicit Float64 setting.

### Spherical grids

HEALPix and cubed-sphere grids support quadrilateral Level-2 footprints on `cpu`, `cuda`, and
`metal`, including standard deviations and centered rolling windows.

```julia
spec = HEALPixGridSpec(nside=2048)
grid_l2(config, spec, time_spec; backend=resolve_backend("cuda"),
        geometry_precision=Float32, compute_std=true, outfile="healpix.nc")
```

```bash
julia --project=. bin/grid.jl l2 --config examples/tropomi_sif.toml \
    --gridType healpix --nside 2048 --backend cuda --geometry-precision float32 \
    --compSTD -o healpix.nc
```

For the cube, pass `CubedSphereGridSpec(Nc=360)` in Julia or use
`--gridType cs --Nc 360` on the CLI. See [cubed-sphere precision](cubed_sphere.md#backends-and-precision)
for measured cell-assignment differences.

Spherical geometry defaults to `Float64` on CUDA and KA CPU. To trade a small
number of cell-assignment differences near edges for potentially faster
execution, pass `geometry_precision=Float32` or `--geometry-precision float32`.
This setting controls the adaptive sample count, unit-vector interpolation,
and cell location. The host also quantizes the corner vectors to Float32
before upload and the kernel orders them in Float32. This geometry setting is
independent of accumulator precision. Metal requires Float32 geometry.

To archive Float64 results exactly, accumulate in Float64 and set
`[output] precision = "float64"`, `significant_digits_mean = "none"`, and
`significant_digits_std = "none"`. On the CLI, KA CPU and CUDA use Float64
accumulators by default. Sequential rolling runs also accumulate in Float64. For tiled
sequential runs in Julia, use Float64 bounds and spacings in `GridSpec` or
`T=Float64` in `HEALPixGridSpec` and `CubedSphereGridSpec`.
Use this configuration for HEALPix 4:1 coarsening or regridding. At fine
`nside`, quantization of nearly equal child means can dominate a small parent
standard deviation.

At `nside=2048`, Float32 shifts HEALPix cell edges by at most about 0.33 m.
For representative OCO and TROPOMI footprints, the combined corner,
interpolation, and locator displacement is about 0.2 m RMS and at most
about 1.1 m, below the typical 50–500 m L2 geolocation uncertainty.
A flipped subpoint can still change occupancy or
shift a mean by an amount comparable to the scene standard deviation in
sparsely grazed perimeter cells (`n ≪ 1`, where `n` is the sum of footprint
weights). In the XCO₂ examples, cells with `n ≳ 2` shifted by about 0.05 ppm
or less.

## Architecture

### Rectangular grids

Level-2 input preparation orders quadrilateral corners on the host. The
rectangular KA index kernels compute each footprint's sample indices and
mark out-of-grid samples. The index kernels can mark single-cell
footprints for a direct update. The quadrilateral kernel can also skip
footprints that are too wide.
Circular sampling counts all in-ellipse samples, including those dropped
at the grid boundary.

Quadrilateral scatter adds one contribution per sub-sample with weight
`1/n²`. Circular scatter counts hits in each cell and adds one contribution
per footprint and cell, dividing the hit count by the in-ellipse sample
count.

### Spherical grids

Spherical gridding filters invalid observations and converts the four corners to
Cartesian unit vectors on the host. The geometry step orders the corners in
the tangent plane, selects the sample count for each footprint, and locates
the bilinear subpoints directly from their normalized vectors. A private
histogram holds eight cells and wide footprints are rescanned in batches of
cell IDs. Each occupied cell receives one contribution per footprint. GPU
kernels fuse geometry and atomic scatter. KA CPU workers count cells in
parallel batches and buffer the counts for scatter in footprint order.

### Accumulation and standard deviation

KA CPU scatters rectangular and spherical contributions serially in
footprint order. GPU backends use atomic updates.

The first pass accumulates weights and weighted values shifted by a fixed
per-variable reference. Each reference is the median of the values finite in
the accumulator precision in the first nonempty batch, or zero if that
variable has no finite values in the batch. Tiled runs choose a reference for
each tile. Rolling runs share the first reference across days.

Float32 runs and runs that request standard deviation record the observed
range of each cell in the first pass. The driver clamps the preliminary
mean to these bounds and rereads each included file. Weighted linear
deviations from that preliminary mean refine it. Float32 runs clamp the
refined mean to the same observed range. When standard deviation
is requested, the replay also accumulates squared deviations, then
subtracts the squared linear residual divided by the weight.
Constant cells with finite moments have zero standard deviation. A
non-finite variance or variance input in accumulator precision produces
fill for the standard deviation. With Float32 accumulators, including on
Metal, the relative accuracy of the mean and standard deviation degrades as
the number of cell contributions grows. Worst-case summation roundoff grows
roughly as `N × eps(Float32)`, where `N` is the number of contributions.
For example, with 100,000 observations, relative mean errors of order
`1e-4` can occur, and standard deviations below that relative scale can be
lost. Float64 accumulators, the default on KA CPU and CUDA, avoid this
Float32 precision limit.
Input files should remain stable for both passes.

Rolling windows cache each day's weight and mean, plus the centered sum
of squared deviations when requested. KA runs that replay do so once per
calendar day and retain a low-order correction to the mean. Float64 stores
the mean relative to the shared reference. Float32 stores the absolute
refined mean so it retains the observed bounds on promotion to Float64.
Centered daily moments merge using Chan's formula in Float64.
Sequential rolling runs also accumulate each day in Float64.

Mean-only Float64 KA runs use a single pass. Their rolling windows add
daily shifted sums in Float64 and restore the shared reference at
finalization. With a distant reference, values with magnitude below about
`eps(Float64) * abs(reference - value)` can round away in the first pass.
The replay handles first-pass drift in Float32, including on Metal.
Float32 precision can still absorb an overlap weight smaller than half an
ulp of the accumulated cell weight. Very small squared deviations can
underflow in Float32.

A new spherical grid implements the [`AbstractSphericalGridSpec`](@ref)
locator for use with both sequential and KA subdivision. Provide angular
cell spacing in radians to enable adaptive sampling.

## Boundary Behavior and Geometry Precision

### Geometry Precision and Sequential Parity

The rectangular KA paths mirror the sequential arithmetic — element types,
subdivision formulas, and rounding steps — so that on KA CPU and CUDA the chosen
oversampling factor `n` and the `(footprint, sample, cell)` assignments match
the sequential path:

- **Fractional indices and automatic `n` (before the kernels)**: Both sequential
  and KA paths compute fractional grid indices directly from the input
  coordinates and `GridSpec{T}` bounds in `promote_type(eltype(coords), T)`, and
  select automatic `n` from those indices before any backend transfer. When both
  the input coordinates and `GridSpec` are `Float32`, `Float32` rounding (such
  as `prevfloat(90f0)` mapping to fractional row `91.0f0` on a `2°` global grid)
  already occurs during this pre-kernel index calculation on both paths.
- **Quadrilateral kernels (KA CPU and CUDA)**: Corner coordinates keep their
  input element type; edge subdivisions are widened to `Float64` baseline points
  and connecting lines are subdivided in `Float64`, then rounded to the grid
  geometry type `T` before `floor`, matching `compute_subpixels!` and
  `floor_indices!`.
- **Circular kernels (KA CPU and CUDA)**: Centers and corners are promoted
  internally to `Float64` to compute symmetric radii and evaluate the sample
  loop, matching `_circle_footprint_axes` and `accumulate_circular_footprint!`.
  Both paths update a cell once per footprint. A single-cell footprint
  contributes exactly one unit of weight.
- **Weights and means**. Both quadrilateral paths update once per sub-sample
  with weight `1/n²`, an exact binary fraction for power-of-two `n`.
  For other values of `n`, the accumulated weights can round differently
  when the accumulator precisions differ. Sequential uses Welford's online
  mean and KA uses reference-shifted sums, so their means can also differ
  at rounding level.
- **Metal exception**: Metal does not support `Float64`. Any `Float64`
  coordinate inputs are converted to `Float32` on transfer, and both
  quadrilateral intermediate baseline interpolation and circular radius/sample
  arithmetic execute in `Float32`. On Metal, footprints or sub-samples right on
  cell boundaries may therefore be assigned to an adjacent cell relative to the
  sequential path.

### Out-of-Grid Sub-Samples and Global Boundaries

- **Drop without renormalization**: On all rectangular paths (sequential and KA,
  quadrilateral and circular), sub-samples that fall outside
  `1..n_lon × 1..n_lat` are dropped without renormalizing the remaining in-grid
  weights (`1/n²` for quadrilaterals, `1/inside_count` for circles). Dropping
  without renormalization keeps regional boundary cells domain-invariant so
  accumulated weights (`"n"`) reflect true fractional footprint overlap and a
  grazing edge footprint cannot bias a boundary cell's mean relative to interior
  footprints.
- **Global edges (`±180°` and poles)**: On global rectangular grids, sub-samples
  beyond `±180°` longitude or beyond `±90°` latitude are dropped rather than
  wrapped across the antimeridian or folded across the poles.
- **Spatial filtering in `grid_l2`**: `apply_filters` requires all four
  footprint corners to lie strictly inside `(lat_min, lat_max) × (lon_min, lon_max)`
  (and longitude span `< 50°`), excluding edge-straddling footprints before
  accumulation.
- **Circular footprint radius**: Circular footprints infer symmetric latitude
  and longitude semi-axes from the farthest corner (`max |corner - center|`
  along each axis), which introduces `≲ 0.1%` edge weight bias for physical
  GOSAT-like viewing geometry.

## Performance Considerations

### Data Transfer

Grid accumulators are allocated on the GPU and stay there across files. Each
file's input data is transferred to the GPU for processing. Input I/O and
host preparation can dominate for workloads with small footprints or coarse
grids.

### Memory

Rectangular GPU memory usage scales with:
- `n_fp × n² × sizeof(Int32) × 2` for index arrays (`ix`, `iy`)
- `n_fp × n_vars × sizeof(A)` for values
- `n_lon × n_lat × (n_vars + 1) × sizeof(A)` for weights and shifted sums
- For replay, two buffers per variable in Float64 standard deviation runs
  and three or four in Float32 runs (mean only or with standard deviation).
  Each buffer takes `n_lon × n_lat × sizeof(A)` bytes. Rectangular
  `grid_l2` groups circular samples in the existing index arrays.

For a typical TROPOMI file (~500k soundings, n=10, 6 variables) in Float64:
- Index arrays: ~400 MB
- Values: ~24 MB
- Grid (360×180): ~3.6 MB for the mean and ~6.2 MB extra for standard deviation

The spherical path uploads `n_fp × 4 × 3` corner components in the geometry
precision and `n_fp × n_vars` values in the accumulator precision per file.
Dense accumulators dominate GPU memory for fine spherical grids. The weights
and first-pass sums use `n_cells × (n_vars + 1) × sizeof(A)` bytes,
with `n_cells = 12 × nside²` for HEALPix and `6 × Nc²` for the cube.
Float64 replay reuses the two
observed-bound buffers for deviations. Float32 replay retains the two
bound buffers and uses a separate buffer for linear deviations and another
for squared deviations when computing standard deviation. At `nside=2048`
with six variables, peak dense Float64 allocations are about 2.82 GB for
mean only or 7.65 GB (7.13 GiB) with standard deviation. Float32 needs about
5.03 GB for mean only or 6.24 GB with standard deviation. Per-file inputs
are additional.
In rolling runs with standard deviation, each cached day at that resolution
uses about 7.65 GB in Float64 or 3.83 GB in Float32. A cached Float32
mean-only day stores two parts of the refined mean plus the weight
(about 2.62 GB). The corresponding Float64 window buffer is about 5.23 GB
for refined means or 7.65 GB with standard deviation.

### Batch Size

`grid_l2` uploads one whole input file at a time. For manual batching of
rectangular quadrilateral footprints, `accumulate_batch!` processes one
batch per call. For fine HEALPix grids, dense accumulators can dominate
GPU memory even with small input files.

## Backend Comparison

| Aspect | Sequential | KA CPU | KA CUDA | KA Metal |
|:-------|:-----------|:-------|:--------|:---------|
| Algorithm | Welford (incremental mean) | Shifted sums | Shifted sums | Shifted sums |
| Parallelism | None | Parallel geometry | GPU-parallel | GPU-parallel |
| STD support | Single-pass | Two-pass | Two-pass | Two-pass |
| Atomic ops | None | None | Backend atomics | Backend atomics |
| Typical hardware | Any CPU | Any CPU | NVIDIA GPU | Apple GPU |
