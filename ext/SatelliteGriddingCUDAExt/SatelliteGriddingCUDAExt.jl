module SatelliteGriddingCUDAExt

using SatelliteGridding
using CUDA
using KernelAbstractions

# When CUDA is loaded, the KernelAbstractions kernels in SatelliteGridding
# automatically work with CuArrays via the CUDABackend().
# This extension provides convenience constructors and GPU-specific optimizations.

"""
    cuda_backend()

Return the CUDA backend for use with `grid_l2` or the dispatch-based `grid`
API. Requires `using CUDA`.

# Example
```julia
using SatelliteGridding, CUDA
grid_l2(config, grid_spec, time_spec; backend=cuda_backend())
```
"""
cuda_backend() = CUDABackend()

end # module
