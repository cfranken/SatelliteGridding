using Test
using SatelliteGridding
using Dates
using Statistics
using NCDatasets
using OrderedCollections
using KernelAbstractions
using Zarr   # loads the SatelliteGriddingZarrExt extension for the Zarr daily-emit tests

@testset "SatelliteGridding.jl" begin
    include("test_oversampling.jl")
    include("test_averaging.jl")
    include("test_ka_kernels.jl")
    include("test_backends.jl")
    include("test_filters.jl")
    include("test_config.jl")
    include("test_ncio.jl")
    include("test_vegetation.jl")
    include("test_synthetic_geometry.jl")
    include("test_cubed_sphere.jl")
    include("test_integration.jl")
    include("test_rolling.jl")
    include("test_zarr_ordering.jl")
end
