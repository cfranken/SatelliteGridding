module SatelliteGriddingMetalExt

using SatelliteGridding
using Metal
using KernelAbstractions

"""
    metal_backend()

Return the Metal backend for use with `grid_l2` or the dispatch-based `grid`
API. Requires `using Metal`.
"""
metal_backend() = MetalBackend()

end # module
