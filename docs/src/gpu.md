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

GPU backends support quadrilateral footprints (`--footprint quad`) and circular
footprints (`--footprint circle`, `CircularFootprintGridding`). The circular path
uses a separate index kernel that masks samples outside the inferred circle or
ellipse before scatter accumulation.

## Architecture

The GPU pipeline consists of three KernelAbstractions kernels that run entirely on
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
weights to compute the mean.

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
file's input data is transferred to the GPU for processing. The bottleneck is
typically I/O (reading NetCDF files) rather than computation.

### Memory

GPU memory usage scales with:
- `n_fp × n² × sizeof(Int32) × 2` for index arrays (`ix`, `iy`)
- `n_fp × n_vars × sizeof(Float32)` for values
- `n_lon × n_lat × (n_vars + 1) × sizeof(Float32)` for grid accumulators

For a typical TROPOMI file (~500k soundings, n=10, 6 variables):
- Index arrays: ~400 MB
- Values: ~12 MB
- Grid (360×180): ~1.5 MB

### Batch Size

Large files may need to be processed in batches to fit in GPU memory. The
`accumulate_batch!` function handles one batch at a time.

## Backend Comparison

| Aspect | Sequential | KA CPU | KA CUDA | KA Metal |
|:-------|:-----------|:-------|:--------|:---------|
| Algorithm | Welford (incremental mean) | Sum-based | Sum-based | Sum-based |
| Parallelism | None | Multi-threaded | GPU-parallel | GPU-parallel |
| STD support | Single-pass | Two-pass | Two-pass | Two-pass |
| Atomic ops | Not needed | Not needed (sequential scatter) | Backend atomics | Backend atomics |
| Typical hardware | Any CPU | Any CPU | NVIDIA GPU | Apple GPU |

The KA CPU backend uses sequential scatter (no atomics) for the accumulation
step, avoiding the overhead of software atomics on CPU. The sorting and subpixel
kernels still run in parallel via KernelAbstractions threading.
