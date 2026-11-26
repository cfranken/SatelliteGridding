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
ellipse before scatter accumulation.

### HEALPix

HEALPix supports quadrilateral Level-2 footprints on `cpu`, `cuda`, and
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

HEALPix geometry defaults to `Float64` on CUDA and KA CPU. To trade a small
number of cell-assignment differences near edges for potentially faster
execution, pass `geometry_precision=Float32` or `--geometry-precision float32`.
This setting controls the adaptive sample count, unit-vector interpolation,
and cell location. The host also quantizes the corner vectors to Float32
before upload and the kernel orders them in Float32. Set
`accumulator_precision=Float32` or `--accumulator-precision float32` to
independently select the precision for weights, shifted sums, deviations,
and rolling moments. CUDA and KA CPU default to Float64 for both settings.
Metal requires Float32 for both and rejects an explicit Float64 setting.

To archive Float64 results exactly, accumulate in Float64 and set
`[output] precision = "float64"`, `significant_digits_mean = "none"`, and
`significant_digits_std = "none"`. On the CLI, HEALPix with `--backend cpu`
or `--backend cuda` uses Float64 accumulators by default. Other CLI modes
accumulate in Float32. For sequential HEALPix in Julia, construct the grid
with `HEALPixGridSpec(; nside=..., T=Float64)`. For other grids, use Float64
bounds and spacings in `GridSpec` or `T=Float64` in `CubedSphereGridSpec`.
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

The rectangular GPU pipeline consists of three KernelAbstractions kernels that run entirely on
the GPU. Quadrilateral and circular footprints share the same scatter pattern,
but use different index-generation kernels.

### Kernel 1: Corner Sorting

```
sort_corners_ccw_ka!(backend, lat_corners, lon_corners)
```

Sorts all footprint corners into CCW order in parallel using a 5-comparator
sorting network. Each thread handles one footprint.

### Kernel 2: Subpixel Index Computation

```
compute_footprint_indices_ka!(backend, ix, iy, skip_flag,
                               lat_corners, lon_corners, n)
```

Computes n×n subpixel grid cell indices per footprint. Each thread handles one
footprint and produces n² index pairs. The `skip_flag` indicates:
- `0` = oversample (n×n subpixels computed)
- `1` = fast path (all corners in same cell, single index)
- `2` = skip (footprint too wide for meaningful oversampling)

For circular footprints:

```
compute_circular_footprint_indices_ka!(backend, ix, iy, inside_count,
                                       skip_flag, center_lat, center_lon,
                                       lat_corners, lon_corners, n)
```

This kernel samples the footprint bounding box, keeps only samples inside the
inferred circle/ellipse, and records `inside_count` so scatter weights are
normalized as `1 / inside_count`.

### Kernel 3: Scatter-Accumulate

```
scatter_accumulate_ka!(backend, grid_sum, grid_weights,
                        ix, iy, skip_flag, values, n, n_vars)
```

Atomically adds weighted values to the grid cells using `@atomic`. On GPU
backends this maps to backend-supported atomic operations. Each thread handles
one footprint and scatters its n² subpixel contributions.

Circular footprints use `scatter_accumulate_circular_ka!`, which skips masked
samples and uses per-footprint normalization.

After processing all files for a time step, `finalize_mean!` divides sums by
weights to compute the mean. Rectangular KA accumulates unshifted `Σw·x`, so
Float32 means can drift for many near-constant, high-offset values.

### Spherical grids

HEALPix filters invalid observations and converts the four corners to
Cartesian unit vectors on the host. A fused kernel orders the corners in the
tangent plane, selects the sample count for each footprint, and locates the
bilinear subpoints directly from their normalized vectors. The kernel holds
a private histogram of eight cells and rescans wide footprints in batches
of cell IDs when necessary. Each occupied cell receives one atomic
contribution per footprint. The private histograms replace the rectangular
path's global subpoint index tables.

The first pass accumulates weights and weighted values shifted by a fixed
per-variable reference. Each reference is the median of the values finite in
the accumulator precision in the first nonempty batch, or zero if that
variable has no finite values in the batch. Tiled runs choose a reference for
each tile. Rolling runs share the first reference across days.

For standard deviation, the first pass bounds each cell's preliminary mean
by the observed range. The driver rereads each included file and accumulates
weighted linear and squared deviations from that cell mean. The linear
residual refines the mean and corrects the sum of squared deviations.
Constant cells with finite moments have zero standard deviation. A
non-finite variance or variance input in accumulator precision produces
fill for the standard deviation. With Float32 accumulators, including on
Metal, a standard deviation below roughly one ulp of the cell mean
(`σ/|mean| ≲ 1e-7`) can be reported as zero. Use Float64 accumulators
(the CUDA and KA CPU default) when such small relative variability matters.
Input files should remain stable for both passes.

Rolling windows run both passes once per calendar day and cache the weight,
reference-shifted mean with its low-order correction, and central sum of
squared deviations. The window combines those moments using differences
between the shifted means and restores the reference at finalization. Sequential
rolling windows use unshifted raw moments. For high-offset tracers such as
XCO₂ with small variability, use `HEALPixGridSpec(nside=...; T=Float64)`
with the sequential backend.

Mean-only runs skip the refinement pass. In either accumulator precision
`A`, values with magnitude below about `eps(A) * abs(reference - value)`
can round away when the distant per-variable reference is subtracted.
`eps(Float64)` is about `2e-16`, and `eps(Float32)` is about `1e-7`.
For example, with a reference of `1f0`, both `1f-8` and `2f-8` can give a
zero mean together in Float32. Enable `compute_std` to recover those
values with the refinement pass. Repeated Float32 additions can also drift
even when differences are representable. For 10,000 equal-weight
`410.1f0` observations with a reference of `380f0`, KA CPU measured a
mean about 0.003 ppm low in Float32. Float32 can also absorb an overlap
weight smaller than half an ulp of the accumulated cell weight. With
standard deviations enabled, the second pass refines the mean on Metal,
subject to Float32 limits on the accumulated weights and moments. Very
small squared deviations can underflow in Float32. Float64 extends the
range and precision of the same arithmetic.

A new spherical grid opts in by defining
`SatelliteGridding._spherical_subdivision(spec) = Val(:unit_vectors)` so
the sequential path uses the same subdivision, and by defining
`SatelliteGridding._device_locator(spec)` to return an isbits value and
`SatelliteGridding.locate(locator, x, y, z)` to return the one-based linear
accumulator index of a normalized vector. The locator must use device-safe
operations for both geometry precisions the grid supports. Implement
`SatelliteGridding._device_cell_spacing(locator, F)` to enable adaptive
subdivision with angular cell spacing in radians. The default uses four
samples per side. The cubed sphere uses the sequential CPU backend.

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
  accumulator type `T` before `floor`, matching `compute_subpixels!` and
  `floor_indices!`.
- **Circular kernels (KA CPU and CUDA)**: Centers and corners are promoted
  internally to `Float64` to compute symmetric radii and evaluate the sample
  loop, matching `_circle_footprint_axes` and `accumulate_circular_footprint!`.
  Single-cell circular footprints take a fast path that adds weight `1` once per
  footprint, whereas the sequential path repeatedly adds the rounded per-sample
  weight (`T(1 / inside_count)`) once per inside sample. Cell assignments remain
  identical, while the sequential weight sum drifts with cell occupancy (for
  example, 1,000 single-cell `Float32` footprints at `n = 20` accumulate
  `1001.98` sequentially versus `1000` on KA). Remaining circular differences
  stay within per-cell tolerance.
- **Accumulators and means**: Grid accumulators (`grid_sum`, `grid_weights`)
  remain in the `GridSpec` element type `T`. Finalized means can still differ at
  floating-point rounding level because the sequential path uses Welford's
  online update whereas KA accumulates weighted sums.
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
- `n_fp × n_vars × sizeof(Float32)` for values
- `n_lon × n_lat × (n_vars + 1) × sizeof(Float32)` for grid accumulators

For a typical TROPOMI file (~500k soundings, n=10, 6 variables):
- Index arrays: ~400 MB
- Values: ~12 MB
- Grid (360×180): ~1.5 MB

The HEALPix path uploads `n_fp × 4 × 3` corner components in the geometry
precision and `n_fp × n_vars` values in the accumulator precision per file.
Dense accumulators dominate GPU memory for fine HEALPix grids. The weights
and first-pass sums use `12 × nside² × (n_vars + 1) × sizeof(A)` bytes,
where `A` is the accumulator precision. Standard deviation adds two
per-variable buffers for the squared and linear deviations. Those buffers
hold the observed bounds in the first pass. The peak dense device
allocation is `12 × nside² × (3 × n_vars + 1) × sizeof(A)` bytes with std.
At `nside=2048` with six variables, mean-only needs about 2.82 GB in
Float64 and std needs about 7.65 GB (7.13 GiB). Float32 needs about
1.41 GB and 3.83 GB respectively. Per-file inputs are additional.
Rolling std also needs about 7.65 GB of Float64 host memory for each cached
day with data and another 7.65 GB for the reusable window workspace at that
resolution and variable count.

### Batch Size

`grid_l2` uploads one whole input file at a time. For manual batching of
rectangular quadrilateral footprints, `accumulate_batch!` processes one
batch per call. For fine HEALPix grids, dense accumulators can dominate
GPU memory even with small input files.

## Backend Comparison

| Aspect | Sequential | KA CPU | KA CUDA | KA Metal |
|:-------|:-----------|:-------|:--------|:---------|
| Algorithm | Welford (incremental mean) | Sum-based | Sum-based | Sum-based |
| Parallelism | None | Multi-threaded | GPU-parallel | GPU-parallel |
| STD support | Single-pass | HEALPix only, two-pass | HEALPix only, two-pass | HEALPix only, two-pass |
| Atomic ops | None | Spherical scatter | Backend atomics | Backend atomics |
| Typical hardware | Any CPU | Any CPU | NVIDIA GPU | Apple GPU |

For rectangular grids, KA CPU uses sequential scatter to avoid the overhead
of software atomics on CPU. Its sorting and subpixel kernels run in parallel
via KernelAbstractions threading. Rectangular KA backends do not support
`compute_std=true`. For HEALPix, KA CPU uses the same fused,
atomic kernel as the GPU backends.
