# GPU Acceleration

`SatelliteGridding.jl` supports GPU acceleration via
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) and
[CUDA.jl](https://github.com/JuliaGPU/CUDA.jl). The same kernel code runs on both
CPU and GPU backends.

## Setup

CUDA.jl is a weak dependency — it's only loaded when you explicitly import it:

```julia
using CUDA
using KernelAbstractions
using SatelliteGridding
```

Verify your GPU is detected:
```julia
CUDA.devices()   # Should list your GPU(s)
```

## Usage

### Julia API

```julia
grid_l2(config, grid_spec, time_spec;
        backend = CUDABackend(),
        outfile = "output.nc")
```

### CLI

```bash
julia --project=. bin/grid.jl l2 \
    --config examples/tropomi_sif.toml \
    --backend cuda \
    -o output.nc
```

## Architecture

The GPU pipeline consists of three KernelAbstractions kernels that run entirely on
the GPU:

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

### Kernel 3: Scatter-Accumulate

```
scatter_accumulate_ka!(backend, grid_sum, grid_weights,
                        ix, iy, skip_flag, values, n, n_vars)
```

Atomically adds weighted values to the grid cells using `@atomic`. On CUDA,
this uses hardware atomic operations. Each thread handles one footprint and
scatters its n² subpixel contributions.

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

| Aspect | Sequential | KA CPU | KA CUDA |
|:-------|:-----------|:-------|:--------|
| Algorithm | Welford (incremental mean) | Sum-based | Sum-based |
| Parallelism | None | Multi-threaded | GPU-parallel |
| STD support | Single-pass | Two-pass | Two-pass |
| Atomic ops | Not needed | Not needed (sequential scatter) | Hardware atomics |
| Typical speedup | 1× | 1–2× | 5–10× |

The KA CPU backend uses sequential scatter (no atomics) for the accumulation
step, avoiding the overhead of software atomics on CPU. The sorting and subpixel
kernels still run in parallel via KernelAbstractions threading.
