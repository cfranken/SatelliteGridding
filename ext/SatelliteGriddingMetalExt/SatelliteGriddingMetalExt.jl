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

SatelliteGridding._default_ka_precision(::MetalBackend) = Float32
SatelliteGridding._validate_ka_precision(::MetalBackend, F, A) =
    SatelliteGridding._validate_metal_precision(F, A)

end # module
