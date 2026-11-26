#!/usr/bin/env julia
# Usage: julia --project=. bench/healpix_gpu.jl
# CUDA_DEVICE=1 selects the second GPU by default. Check nvidia-smi first.
# Set SKIP_CUDA=1 to measure the sequential baseline alone.

using SatelliteGridding
using KernelAbstractions
using Random
using Printf

const SG = SatelliteGridding
const EARTH_RADIUS_KM = 6371.0
const ORBIT_SEGMENT_LINES = 2000
const WARMUP_FOOTPRINTS = 1024
const SKIP_CUDA = get(ENV, "SKIP_CUDA", "0") == "1"

if !SKIP_CUDA
    using CUDA
end

"""
Generate consecutive swath rows on short great-circle orbit segments.
Offsets in the along-track and cross-track directions are in kilometres.
"""
function synthetic_day(rng, footprints, batch, cross_track, along_km, nadir_km, edge_km)
    batches = @NamedTuple{lat_bounds::Matrix{Float64}, lon_bounds::Matrix{Float64},
                          values::Matrix{Float32}}[]
    for first in 1:batch:footprints
        count = min(batch, footprints - first + 1)
        lat, lon = Matrix{Float64}(undef, count, 4), Matrix{Float64}(undef, count, 4)
        values = rand(rng, Float32, count, 1)

        # Integrate the across-track width profile. A fourth-power increase
        # puts the broadest footprints near the edges and gives the synthetic
        # TROPOMI swath a width of approximately 2600 km.
        detector = range(-1.0, 1.0; length=cross_track + 1)
        edges = (cross_track / 2) .* (nadir_km .* detector .+
                                      ((edge_km - nadir_km) / 5) .* detector.^5)
        origin = heading = across = (0.0, 0.0, 0.0)
        @inbounds for row in 0:cld(count, cross_track)-1
            offset = row % ORBIT_SEGMENT_LINES
            if iszero(offset)
                lambda = 2pi * rand(rng)
                phi = asin(2rand(rng) - 1)
                sin_lon, cos_lon = sincos(lambda)
                sin_lat, cos_lat = sincos(phi)
                origin = (cos_lat * cos_lon, cos_lat * sin_lon, sin_lat)
                north = (-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat)
                east = (-sin_lon, cos_lon, 0.0)
                inclination = rand(rng, (-1, 1)) * deg2rad(10 + 10rand(rng))
                heading = cos(inclination) .* north .+ sin(inclination) .* east
                across = SG._sphere_cross(heading, origin)
            end
            for (corner, (a, b)) in enumerate(((0, 0), (0, 1), (1, 1), (1, 0)))
                along = ((offset + a) * along_km) / EARTH_RADIUS_KM
                along_position = cos(along) .* origin .+ sin(along) .* heading
                for col in 1:min(cross_track, count - row * cross_track)
                    k = row * cross_track + col
                    sweep = edges[col + b] / EARTH_RADIUS_KM
                    x, y, z = cos(sweep) .* along_position .+ sin(sweep) .* across
                    lat[k, corner] = atand(z, hypot(x, y))
                    lon[k, corner] = atand(y, x)
                end
            end
        end
        push!(batches, (lat_bounds=lat, lon_bounds=lon, values=values))
    end
    return batches
end

function time_sequential(spec, batches)
    shape = accumulator_shape(spec)
    data, weights = zeros(Float32, shape..., 1), zeros(Float32, shape)
    n = SG._HEALPIX_MAX_OVERSAMPLE
    cells = Vector{Tuple{Int,Int}}(undef, n^2)
    seconds = @elapsed for observations in batches
        SG._accumulate_unit_vector_footprints!(data, nothing, weights, false, spec,
            observations.lat_bounds, observations.lon_bounds, observations.values, nothing, cells)
    end
    return seconds, sum(Float64, weights)
end

function time_cuda(backend, spec, batches; geometry_precision, accumulator_precision)
    options = SG._spherical_ka_options(spec, backend, geometry_precision, accumulator_precision)
    F, A = options.geometry_precision, options.accumulator_precision
    shape = accumulator_shape(spec)
    data = KernelAbstractions.zeros(backend, A, shape..., 1)
    weights = KernelAbstractions.zeros(backend, A, shape)
    reference = SG._SphericalKAReference(backend, A, size(data, 3))
    KernelAbstractions.synchronize(backend)
    seconds = @elapsed for observations in batches
        SG._accumulate_spherical_batch!(backend, data, weights, nothing, options.locator, F,
            observations.lat_bounds, observations.lon_bounds, observations.values, nothing, reference)
    end
    return seconds, sum(Float64, Array(weights))
end

function main()
    backend = if SKIP_CUDA
        nothing
    else
        CUDA.functional(true)
        device = parse(Int, get(ENV, "CUDA_DEVICE", "1"))
        CUDA.device!(device)
        println("CUDA device $device: ", CUDA.name(CUDA.device()))
        resolve_backend("cuda")
    end
    println("Accumulation only, including host corner preparation, input upload, and synchronization.")
    println("One variable and adaptive n. The sequential baseline uses Float32 accumulators.")
    println("Data generation, dense allocation, finalization, output download, and NetCDF I/O are outside timing.")
    cuda_precisions = ((nothing, nothing), (Float32, nothing), (nothing, Float32), (Float32, Float32))
    rng = MersenneTwister(413)
    warmup = synthetic_day(rng, WARMUP_FOOTPRINTS, WARMUP_FOOTPRINTS, 450, 5.5, 3.5, 14.5)
    spec = HEALPixGridSpec(nside=64)
    time_sequential(spec, warmup)
    if backend !== nothing
        for (geometry_precision, accumulator_precision) in cuda_precisions
            time_cuda(backend, spec, warmup; geometry_precision, accumulator_precision)
        end
    end
    for (product, footprints, batch, detectors, along, nadir, edge) in (
        ("TROPOMI-like", 3_000_000, 500_000, 450, 5.5, 3.5, 14.5),
        ("OCO-like", 1_000_000, 250_000, 8, 2.25, 1.3, 1.3),
    )
        batches = synthetic_day(rng, footprints, batch, detectors, along, nadir, edge)
        for nside in (64, 2048)
            spec = HEALPixGridSpec(; nside)
            cpu_seconds, weight = time_sequential(spec, batches)
            @assert isapprox(weight, footprints; rtol=1e-5)
            @printf("%-13s %9d footprints, nside=%4d | sequential: %8.3f s\n",
                    product, footprints, nside, cpu_seconds)
            if backend !== nothing
                for (geometry_precision, accumulator_precision) in cuda_precisions
                    seconds, weight = time_cuda(backend, spec, batches;
                                                geometry_precision, accumulator_precision)
                    @assert isapprox(weight, footprints; rtol=1e-5)
                    options = SG._spherical_ka_options(spec, backend, geometry_precision, accumulator_precision)
                    @printf("  CUDA geometry=%s accumulators=%s%s | %8.3f s  (%5.2f× speedup)\n",
                            string(options.geometry_precision), string(options.accumulator_precision),
                            accumulator_precision === nothing ? " (default)" : "", seconds, cpu_seconds / seconds)
                end
            end
            flush(stdout)
        end
    end
    return nothing
end

main()
