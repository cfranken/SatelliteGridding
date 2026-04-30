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
