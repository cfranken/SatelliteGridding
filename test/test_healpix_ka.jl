using Random
using SHA

# Seam-crossing, broad, polar, pole-enclosing, and face/cell-straddling footprints.
function _healpix_ka_footprints(::Type{T}) where {T}
    geo_lat = [-2.2 -2.1 1.9 2.6; -1.4 -1.3 2.0 2.2; 3.0 3.1 7.6 7.0;
           4.0 4.0 8.0 8.0; 86.0 87.0 89.0 89.2; -0.02 -0.02 0.02 0.02;
           89.8 89.8 89.8 89.8; 38.1 38.1 39.8 39.8; 0.05 0.05 0.2 0.2]
    geo_lon = [32.1 49.3 47.7 34.0; 33.2 50.0 46.8 37.0; 176.1 -176.2 -177.0 177.4;
           177.2 -178.0 -177.3 177.4; 127.0 -137.0 -101.0 113.0;
           -10.02 -9.98 -9.98 -10.02; 45 135 225 315; -0.2 0.2 0.2 -0.2; 0.02 0.5 0.5 0.02]
    vals = T[1.1 9.2; -4.3 6.4; 13.5 -2.6; 31.7 5.8; 0.9 4.1; 8.2 1.3;
             5.0 2.0; 1.4 3.0; 7.2 2.8]
    return geo_lat, geo_lon, vals
end

function _healpix_accumulator_hash(T, nside, n)
    lat, lon, values = _healpix_ka_footprints(T)
    spec = HEALPixGridSpec(; nside, T)
    shape = accumulator_shape(spec)
    data, m2, weights = zeros(T, shape..., 2), zeros(T, shape..., 2), zeros(T, shape)
    observations = (lat_bounds=lat, lon_bounds=lon, values=values)
    SatelliteGridding._accumulate_l2_spherical!(Val(:unit_vectors),
        data, m2, weights, true, spec, observations, n)
    bits = reinterpret(T == Float32 ? UInt32 : UInt64, vcat(vec(data), vec(m2), vec(weights)))
    return bytes2hex(sha256(reinterpret(UInt8, htol.(bits))))
end

function _write_healpix_ka_inputs(dir)
    lat, lon, values = _healpix_ka_footprints(Float64)
    # Shift the deliberately edge-aligned hash fixture slightly for the
    # precision comparison, retaining the symmetric pole-enclosing footprint.
    lat[[1:6; 8:9], :] .+= 0.007
    lon[[1:6; 8:9], :] .+= 0.013
    # Each input has a rejected non-finite value and a rejected latitude.
    lat = vcat(lat, lat[1:2, :])
    lon = vcat(lon, lon[1:2, :])
    lat[end, 1] = 91
    values = vcat(values, [NaN 1.0; 123.0 456.0])
    for (day, orbit, offset) in ((1, 1, 0), (1, 2, 1), (2, 1, 3), (3, 1, -2))
        path = joinpath(dir, "ka_2020-01-0$(day)_$orbit.nc")
        Dataset(path, "c") do ds
            defDim(ds, "sounding", size(lat, 1))
            defDim(ds, "corner", 4)
            defVar(ds, "lat", Float64, ("sounding",))[:] = lat[:, 1]
            defVar(ds, "lon", Float64, ("sounding",))[:] = lon[:, 1]
            defVar(ds, "lat_bounds", Float64, ("sounding", "corner"))[:, :] = lat
            defVar(ds, "lon_bounds", Float64, ("sounding", "corner"))[:, :] = lon
            defVar(ds, "signal", Float64, ("sounding",))[:] = values[:, 1] .+ offset
            defVar(ds, "other", Float64, ("sounding",))[:] = values[:, 2] .- offset
        end
    end
    return DataSourceConfig(Dict("lat" => "lat", "lon" => "lon",
              "lat_bnd" => "lat_bounds", "lon_bnd" => "lon_bounds"),
        OrderedDict("signal" => "signal", "other" => "other"),
        FilterRule[], "ka_YYYY-MM-DD_*.nc", dir)
end

function _compare_healpix_ka_files(reference, actual, expected_weights, std)
    Dataset(reference) do cpu
        Dataset(actual) do ka
            w_cpu, w_ka = cpu["n"][:, :], ka["n"][:, :]
            hits = w_cpu .> 0
            @test hits == (w_ka .> 0)
            @test all(isapprox.(w_cpu[hits], w_ka[hits]; rtol=2e-6))
            @test vec(sum(w_ka; dims=1)) ≈ expected_weights rtol=2e-6
            @test cpu["window_start"][:] == ka["window_start"][:]
            @test cpu["window_end"][:] == ka["window_end"][:]
            for name in ("signal", "other")
                @test all(isapprox.(cpu[name][:, :][hits], ka[name][:, :][hits];
                                    rtol=2e-6, atol=1e-6))
                std_name = name * "_std"
                @test haskey(ka, std_name) == std
                if std
                    @test all(isapprox.(cpu[std_name][:, :][hits], ka[std_name][:, :][hits];
                                        rtol=2e-6, atol=1e-6))
                end
            end
        end
    end
end

function _test_healpix_ka_end_to_end(backend; compare_float32=false)
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        config = _write_healpix_ka_inputs(dir)
        d1, d2 = DateTime(2020, 1, 1), DateTime(2020, 1, 2)
        tiled = TimeSpec(d1, d2, Day(1))
        rolling = TimeSpec(d1, d2, Day(1); sample_step=Day(1), window_halfwidth_days=1)
        seq_path = joinpath(dir, "seq.nc")
        ka_path = joinpath(dir, "ka.nc")
        for (nside, T, times, n, std, weights) in (
            (8, Float32, tiled, 8, false, [18, 9]),
            (8, Float64, tiled, 9, true, [18, 9]),
            (128, Float64, tiled, nothing, true, [18, 9]),
            (8, Float64, rolling, nothing, true, [27, 36]),
            (8, Float32, rolling, 8, false, [27, 36]),
        )
            spec = HEALPixGridSpec(; nside, T)
            opts = (; n_oversample=n, compute_std=std, progress=false)
            grid_l2(config, spec, times; outfile=seq_path, opts...)
            grid_l2(config, spec, times; backend, outfile=ka_path, opts...)
            _compare_healpix_ka_files(seq_path, ka_path, weights, std)
            if compare_float32
                grid_l2(config, spec, times; backend, geometry_precision=Float32,
                        outfile=ka_path, opts...)
                _compare_healpix_ka_files(seq_path, ka_path, weights, std)
            end
        end
    end
end

function _healpix_ka_signal_config(dir, prefix)
    return DataSourceConfig(
        Dict("lat" => "lat", "lon" => "lon", "lat_bnd" => "lat_bounds", "lon_bnd" => "lon_bounds"),
        OrderedDict("signal" => "signal"), FilterRule[], "$(prefix)_YYYY-MM-DD.nc", dir)
end

function _write_healpix_ka_signal_input(path, rectangles, values)
    Dataset(path, "c") do ds
        defDim(ds, "sounding", length(values))
        defDim(ds, "corner", 4)
        defVar(ds, "signal", eltype(values), ("sounding",))[:] = values
        for (name, lo, hi) in (("lon", first.(rectangles), getindex.(rectangles, 2)),
                               ("lat", getindex.(rectangles, 3), last.(rectangles)))
            corners = name == "lat" ? hcat(lo, lo, hi, hi) : hcat(lo, hi, hi, lo)
            defVar(ds, name, Float64, ("sounding",))[:] = (lo .+ hi) ./ 2
            defVar(ds, name * "_bounds", Float64, ("sounding", "corner"))[:, :] = corners
        end
    end
end

function _test_healpix_moment_precision(backend, accumulator_precision)
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        spec = HEALPixGridSpec(nside=8, T=Float32)
        lon, lat = pix2ang_nest(8, 100)
        rect = (lon - 0.01, lon + 0.01, lat - 0.01, lat + 0.01)
        grazes = (-0.03, 0.04, -0.025, 0.02)
        date = DateTime(2020, 1, 1)
        tiled = TimeSpec(date, date, Day(2))
        rolling = TimeSpec(date, date, Day(1); sample_step=Day(1), window_halfwidth_days=1)
        output = joinpath(dir, "output.nc")
        rng = MersenneTwister(510)
        narrow = 410f0 .+ 0.05f0 .* randn(rng, Float32, 10_000)
        for (prefix, values) in (("constant", fill(400.1f0, 4096)), ("narrow", narrow))
            midpoint = length(values) ÷ 2
            _write_healpix_ka_signal_input(joinpath(dir, "$(prefix)_2020-01-01.nc"),
                            fill(rect, midpoint), values[1:midpoint])
            # This separate footprint supplies fractional-weight single-footprint cells.
            _write_healpix_ka_signal_input(joinpath(dir, "$(prefix)_2020-01-02.nc"),
                            [fill(rect, midpoint); [grazes]],
                            [values[midpoint+1:end]; 397.3f0])
            config = _healpix_ka_signal_config(dir, prefix)
            reference = Float64.(values)
            for times in (tiled, rolling)
                grid_l2(config, spec, times; backend, accumulator_precision, compute_std=true,
                        n_oversample=7, outfile=output, progress=false)
                Dataset(output) do ds
                    n = ds["n"][:, 1]
                    signal = ds["signal"][:, 1]
                    std_var = ds["signal_std"]
                    @test n[101] == length(values)
                    @test abs(signal[101] - mean(reference)) <= eps(Float32(mean(reference)))
                    if prefix == "constant"
                        @test std_var[101, 1] <= 1e-5
                    else
                        @test std_var[101, 1] ≈ std(reference; corrected=false) rtol=1e-3
                    end
                    if times === tiled
                        singles = (0 .< n) .& (n .< 1)
                        @test count(singles) > 1
                        @test all(iszero, std_var[:, 1][singles])
                        @test occursin("Zero for a cell touched by a single footprint in tiled mode",
                                       std_var.attrib["comment"])
                    else
                        @test occursin("centered daily moments", std_var.attrib["comment"])
                    end
                end
            end
        end
    end
end

function _healpix_cell_rectangle(nside, pix)
    lon, lat = pix2ang_nest(nside, pix)
    return (lon - 0.01, lon + 0.01, lat - 0.01, lat + 0.01)
end

function _test_healpix_accumulator_cases(backend, accumulator_precision)
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        spec = HEALPixGridSpec(nside=8)
        target = _healpix_cell_rectangle(8, 100)
        other = _healpix_cell_rectangle(8, 500)
        d = DateTime(2020, 1, 1)
        tiled = TimeSpec(d, d, Day(2))
        rolling = TimeSpec(d, d, Day(1); sample_step=Day(1), window_halfwidth_days=1)
        output = joinpath(dir, "output.nc")
        cases = (
            ("prepended", 380f0, fill(400.1f0, 4096), true),
            ("distant", 380f0, fill(410.1f0, 10_000), false),
            ("small", 1f8, Float32[1, 2], false),
            ("tiny", 1f0, Float32[1f-8, 2f-8], false),
            ("sub_ulp_std", 400f0, [fill(445.1f0, 9999); 445.1f0 - 0.003f0], false),
            ("overflow", 1f20, Float32[1f20, 3f20], false),
            ("sum_overflow", 1f18, repeat(Float32[1f18, 3f18], 256), false),
        )
        for (name, anchor, values, same_batch) in cases
            @testset "$name" begin
                first = joinpath(dir, "$(name)_2020-01-01.nc")
                second = joinpath(dir, "$(name)_2020-01-02.nc")
                if same_batch
                    _write_healpix_ka_signal_input(first, vcat([other], fill(target, length(values))),
                                                   [anchor; values])
                else
                    _write_healpix_ka_signal_input(first, [other], [anchor])
                    _write_healpix_ka_signal_input(second, fill(target, length(values)), values)
                end
                config = _healpix_ka_signal_config(dir, name)
                wide = Float64.(values)
                expected_mean = round(Float32(mean(wide)); sigdigits=8)
                # Match the NetCDF writer's six significant digits. For the
                # 512-observation case the mathematical std is at a Float32 tie.
                expected_std = round(Float32(std(wide; corrected=false)); sigdigits=6)
                single = accumulator_precision === Float32
                overflow = name in ("overflow", "sum_overflow")
                for (times, compute_std) in ((tiled, false), (tiled, true), (rolling, true))
                    grid_l2(config, spec, times; backend, accumulator_precision, compute_std,
                            n_oversample=1, outfile=output, progress=false)
                    Dataset(output) do ds
                        @test ds["n"][101, 1] == length(values)
                        found_mean = ds["signal"][101, 1]
                        if !single || compute_std || name == "prepended"
                            if single && overflow
                                @test found_mean ≈ expected_mean rtol=2e-6
                            else
                                @test found_mean === expected_mean
                            end
                        elseif name == "distant"
                            # The fixed 380 ppm reference leaves a 30.1 ppm shifted sum.
                            @test abs(found_mean - expected_mean) < 0.005f0
                        end
                        if compute_std
                            deviation = ds["signal_std"]
                            @test occursin("non-finite variance or variance input", deviation.attrib["comment"])
                            if single && overflow
                                @test ismissing(deviation[101, 1])
                            elseif single && name == "sub_ulp_std"
                                # Float32 replay can cancel variance at this scale.
                                @test deviation[101, 1] === 0f0
                            else
                                @test deviation[101, 1] === expected_std
                            end
                        end
                    end
                end
            end
        end
        if accumulator_precision === Float32
            # A finite Float64 input that overflows on conversion must still
            # contribute weight and give fill for that cell's standard deviation.
            input = joinpath(dir, "input_overflow_2020-01-01.nc")
            _write_healpix_ka_signal_input(input, [other, target, target], [3.0, 1e40, 2.0])
            config = _healpix_ka_signal_config(dir, "input_overflow")
            grid_l2(config, spec, tiled; backend, accumulator_precision, compute_std=true,
                    n_oversample=1, outfile=output, progress=false)
            Dataset(output) do ds
                @test ds["n"][101, 1] == 2f0
                @test ismissing(ds["signal_std"][101, 1])
                @test ds["signal_std"][501, 1] == 0f0
            end
        else
            # The source precision survives until accumulation even though
            # the grid and the serialized payload are Float32.
            input = joinpath(dir, "wide_input_2020-01-01.nc")
            _write_healpix_ka_signal_input(input, [target, target], [1 - 1e-9, 1 + 1e-9])
            config = _healpix_ka_signal_config(dir, "wide_input")
            grid_l2(config, spec, tiled; backend, accumulator_precision, compute_std=true,
                    n_oversample=1, outfile=output, progress=false)
            Dataset(output) do ds
                @test ds["signal"][101, 1] === 1f0
                @test ds["signal_std"][101, 1] === 1f-9
            end
        end
    end
end

function _test_healpix_rolling_offset(backend, accumulator_precision)
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        spec = HEALPixGridSpec(nside=8)
        target = _healpix_cell_rectangle(8, 100)
        other = _healpix_cell_rectangle(8, 500)
        early = DateTime(2019, 1, 1)
        late = early + Day(365)
        times = TimeSpec(early, late, Day(1); sample_step=Day(365), window_halfwidth_days=1)
        config = _healpix_ka_signal_config(dir, "offset")
        first_day = repeat(Float32[444.93, 444.99], 50)
        second_day = repeat(Float32[445.01, 445.07], 50)
        for (date, rect, values) in (
            (early, other, Float32[400]),
            (late - Day(1), target, first_day),
            (late, target, second_day),
        )
            path = joinpath(dir, "offset_$(Date(date)).nc")
            _write_healpix_ka_signal_input(path, fill(rect, length(values)), values)
        end
        output = joinpath(dir, "output.nc")
        grid_l2(config, spec, times; backend, accumulator_precision, compute_std=true,
                n_oversample=1, outfile=output, progress=false)
        Dataset(output) do ds
            @test ds["signal"][501, 1] == 400f0
            @test ds["n"][101, 2] == 200f0
            @test ds["signal"][101, 2] == 445f0
            expected = std(Float64.(vcat(first_day, second_day)); corrected=false)
            @test ds["signal_std"][101, 2] ≈ expected rtol=1e-4
        end
    end
end

function _test_healpix_rolling_mean_correction(backend, accumulator_precision)
    A = something(accumulator_precision, Float64)
    @testset "Rolling sub-ulp daily means ($A)" begin
        scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
        mktempdir(scratch) do dir
            spec = HEALPixGridSpec(nside=8, T=A)
            rect = _healpix_cell_rectangle(8, 100)
            d = DateTime(2020, 1, 1)
            x = A(400)
            u = eps(x)
            for (date, values) in ((d, [x - u, x]), (d + Day(1), [x, x + u]))
                path = joinpath(dir, "correction_$(Date(date)).nc")
                _write_healpix_ka_signal_input(path, fill(rect, length(values)), values)
            end
            config = _healpix_ka_signal_config(dir, "correction")
            output = joinpath(dir, "output.nc")
            tiled = TimeSpec(d, d, Day(2))
            rolling = TimeSpec(d, d, Day(1); sample_step=Day(1), window_halfwidth_days=1)
            tiled_std, rolling_std = map((tiled, rolling)) do times
                grid_l2(config, spec, times; backend, accumulator_precision, compute_std=true,
                        n_oversample=1, outfile=output, progress=false)
                Dataset(output) do ds
                    @test ds["n"][101, 1] == 4
                    @test ds["signal"][101, 1] == x
                    ds["signal_std"][101, 1]
                end
            end
            @test tiled_std == Float32(round(u / sqrt(A(2)); sigdigits=6))
            @test abs(rolling_std - tiled_std) <= eps(tiled_std)
        end
    end
end

function _test_healpix_tiny_weight(backend, accumulator_precision)
    S = SatelliteGridding
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        spec = HEALPixGridSpec(nside=8, T=Float32)
        date = DateTime(2020, 1, 1)
        tiled = TimeSpec(date, date, Day(1))
        rolling = TimeSpec(date, date, Day(1); sample_step=Day(1), window_halfwidth_days=0)
        output = joinpath(dir, "tiny.nc")
        n = 4096
        # Only one sample of this rectangle falls into linear cell 306.
        grazing = (10.0, 10.1, 8.419654254143046, 8.519654254143045)
        lon, lat = pix2ang_nest(8, 100)
        distant = (lon - 0.01, lon + 0.01, lat - 0.01, lat + 0.01)
        input = joinpath(dir, "grazing_2020-01-01.nc")
        config = _healpix_ka_signal_config(dir, "grazing")
        for (rectangles, values) in (([grazing], Float32[12]),
                                     ([distant, grazing], Float32[9, 12]))
            # The second case includes a separate value so the batch median
            # differs from 12. A raw shifted sum is especially small here.
            _write_healpix_ka_signal_input(input, rectangles, values)
            for times in (tiled, rolling)
                # Exercise the original reproduction above the new public cap.
                Base.invokelatest(S._grid_l2, config, spec, times; backend,
                    effective_n_oversample=n, footprint_method=SubpixelGridding(),
                    compute_std=true, outfile=output, progress=false, geometry_precision=nothing,
                    accumulator_precision, zarr_out=nothing, zarr_meta=nothing, zarr_epoch=S.ZARR_EPOCH,
                    keep_going=false, command=nothing)
                Dataset(output) do ds
                    @test ds["n"][306, 1] == inv(Float32(n^2)) < eps(Float32)
                    @test ds["signal"][306, 1] == 12f0
                    @test ds["signal_std"][306, 1] == 0f0
                end
            end
        end
    end
end

@noinline function _healpix_int32_histogram(corners, normal, rescan::Val{R}) where {R}
    return SatelliteGridding._footprint_cell_counts(_HemisphereTestLocator(), corners,
                                                     normal, 32, 0, rescan)
end

function _healpix_int32_histogram_allocations(rescan::Val{R}) where {R}
    S = SatelliteGridding
    corners, normal = S._sphere_ordered_corners(
        map(latlon -> S._sphere_unit_vector(latlon...), ((-5, 5), (-5, 10), (37, 10), (37, 5)))...)
    ids, counts, occupied = _healpix_int32_histogram(corners, normal, rescan)
    @test ids isa NTuple{8,Int}
    @test occupied == 2
    @test sum(counts) == 32^2
    return @allocated _healpix_int32_histogram(corners, normal, rescan)
end

# The diagnostic uses the sequential corner order for both precisions. The
# device still computes its tangent frame and adaptive count in Float32.
@kernel function _healpix_test_sample_ids!(ids, sample_counts, @Const(vertices), locator, spacing, n)
    k = @index(Global, Linear)
    S = SatelliteGridding
    @inbounds corners = (vertices[k, 1], vertices[k, 2], vertices[k, 3], vertices[k, 4])
    ordered, normal = S._sphere_ordered_corners(corners...)
    @inbounds sample_counts[k] = S._unit_vector_samples(ordered, spacing)
    @inbounds for j in 1:n
        left, right = S._sphere_subpoint_edges(corners, j, n)
        for i in 1:n
            x, y, z = S._sphere_subpoint(left, right, normal, i, n)
            ids[i + (j - 1) * n, k] = S.locate(locator, x, y, z)
        end
    end
end

function _healpix_test_precision(backend)
    S = SatelliteGridding
    rng = MersenneTwister(42)
    n, footprints, nside = 8, 2048, 2048
    lon = 360 .* rand(rng, footprints) .- 180
    lat = asind.(2 .* rand(rng, footprints) .- 1)
    vertices = Matrix{NTuple{3,Float32}}(undef, footprints, 4)
    ref = Matrix{Int64}(undef, n^2, footprints)
    ref_n = Vector{Int}(undef, footprints)
    spec = HEALPixGridSpec(; nside)
    locator = S._device_locator(spec)
    spacing = S._healpix_cell_spacing(spec)
    offsets = ((-0.018, -0.041), (-0.018, 0.041), (0.018, 0.041), (0.018, -0.041))
    for k in 1:footprints
        corners, normal = S._sphere_ordered_corners(
            map(((dy, dx),) -> S._sphere_unit_vector(lat[k] + dy, lon[k] + dx / cosd(lat[k])),
                offsets)...)
        vertices[k, :] .= map(v -> map(Float32, v), corners)
        ref_n[k] = S._unit_vector_samples(corners, spacing)
        for j in 1:n
            left, right = S._sphere_subpoint_edges(corners, j, n)
            for i in 1:n
                x, y, z = S._sphere_subpoint(left, right, normal, i, n)
                col, row = S._unit_vector_cell(spec, locator, x, y, z)
                ref[i + (j - 1) * n, k] = col + (row - 1) * nside^2
            end
        end
    end
    inputs = S._spherical_device_input(backend, vertices)
    ids = KernelAbstractions.allocate(backend, Int64, size(ref))
    counts = KernelAbstractions.allocate(backend, Int, size(ref_n))
    _healpix_test_sample_ids!(backend, 64)(ids, counts, inputs, locator,
        S._device_cell_spacing(locator, Float32), n; ndrange=footprints)
    KernelAbstractions.synchronize(backend)
    mismatches = count(pair -> first(pair) != last(pair), zip(ref, Array(ids)))
    fraction = mismatches / length(ref)
    n_changes = count(ref_n .!= Array(counts))
    @info "HEALPix Float32 precision" backend nside mismatches samples=length(ref) fraction n_changes
    @test fraction < 0.001
    @test n_changes <= 1
end

@testset "HEALPix KA backends" begin
    @testset "Staged sequential accumulator bits" begin
        expected = (
            ("e936b5876c28890b0ef50d8a9a7738fc0148aedb9e8f8b65fe9f413f55c84e90",
             "152c5eb90a267e868182b279a2b995cbf29ce269a86e27363272703975b2d51a",
             "bcdab8a568776a67bef7f4c82d8b458daf8d1a294554120ad73b165f099e249a",
             "c3095da1146abd19802da8510577a6fc9b25ffb121d5ad7f61968efed444ba18"),
            ("c26d81b7706983933fc9d9fb5a6859bf888077cbb789222fa24211743cee88f6",
             "17b8aecea53e50c6c2882f98eeb9af888ef7906bbab9a42c0d857dd6941cb6ce",
             "9c08aa99a81d0d32135a42f8cee5b35263830abfffaa5c9df1a02a66741f421c",
             "75bf9cf2241a04a5ba79f082e6a1464ccdd6f99d8377303884a8caafa513bb5f"),
        )
        for (n, hashes) in zip((7, nothing), expected)
            for ((T, side), hash) in zip(Iterators.product((Float32, Float64), (8, 128)),
                                        (hashes[1], hashes[3], hashes[2], hashes[4]))
                @test _healpix_accumulator_hash(T, side, n) == hash
            end
        end
    end

    @testset "Vector locator near Float32 poles" begin
        side = 8192
        for pix in (side^2 - 2, 4side^2 - 3, 8side^2 + 1, 11side^2 + 2)
            lon, lat = pix2ang_nest(side, pix)
            x, y, z = map(Float32, SatelliteGridding._sphere_unit_vector(lat, lon))
            @test abs(z) == 1f0
            @test vec2pix_nest(side, x, y, z) == pix
        end
    end

    @testset "Float32 equatorial edge precision" begin
        S = SatelliteGridding
        side = 2048
        # Position the ascending edge at each longitude, close to the equator.
        for lon in (10.0, 100.0, 200.0, 355.0)
            phase = side * (0.5 + lon / 90)
            lat = asind((phase - round(phase)) / (0.75 * side))
            delta = rad2deg(0.15 / (S._HEALPIX_EARTH_RADIUS * cosd(lat)))
            west, east = ang2pix_nest(side, lon - delta, lat), ang2pix_nest(side, lon + delta, lat)
            @test west != east
            for (sample, expected) in ((lon - delta, west), (lon + delta, east))
                vector = map(Float32, S._sphere_unit_vector(lat, sample))
                @test vec2pix_nest(side, vector...) == expected
            end
        end
    end

    @testset "KA CPU" begin
        backend = resolve_backend("cpu")
        @testset "Largest public oversampling override" begin
            scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
            mktempdir(scratch) do dir
                lon, lat = pix2ang_nest(8, 100)
                _write_healpix_ka_signal_input(joinpath(dir, "max_2020-01-01.nc"),
                    [(lon - 0.01, lon + 0.01, lat - 0.01, lat + 0.01)], Float32[12])
                date = DateTime(2020, 1, 1)
                config = _healpix_ka_signal_config(dir, "max")
                output = joinpath(dir, "output.nc")
                for target in (nothing, backend)
                    grid_l2(config, HEALPixGridSpec(nside=8), TimeSpec(date, date, Day(1));
                            backend=target, n_oversample=1024, outfile=output, progress=false)
                    Dataset(output) do ds
                        @test ds["n"][101, 1] == 1f0
                        @test ds["signal"][101, 1] == 12f0
                    end
                end
            end
        end
        @test _healpix_int32_histogram_allocations(Val(false)) == 0
        @test _healpix_int32_histogram_allocations(Val(true)) == 0
        _test_healpix_ka_end_to_end(backend; compare_float32=true)
        for precision in (nothing, Float32)
            _test_healpix_moment_precision(backend, precision)
            _test_healpix_accumulator_cases(backend, precision)
            _test_healpix_rolling_offset(backend, precision)
            _test_healpix_rolling_mean_correction(backend, precision)
            _test_healpix_tiny_weight(backend, precision)
        end
        _healpix_test_precision(backend)
        @test parse_l2_args(["--config", "input.toml", "--geometry-precision", "float32",
                             "--backend", "cpu"])["geometry-precision"] == "float32"
        @test parse_l2_args(["--config", "input.toml", "--accumulator-precision", "float32",
                             "--backend", "cpu"])["accumulator-precision"] == "float32"
        spec = HEALPixGridSpec(nside=8)
        @test_throws ArgumentError("Metal requires Float32 accumulator_precision") resolve_backend(
            "metal"; grid_spec=spec, accumulator_precision=Float64)
        @test_throws ArgumentError("accumulator_precision must be Float32 or Float64") resolve_backend(
            "cuda"; grid_spec=spec, accumulator_precision=Float16)
        @test_throws ArgumentError("accumulator_precision requires a spherical grid and a KA backend") resolve_backend(
            "sequential"; grid_spec=spec, accumulator_precision=Float32)
        @test_throws ArgumentError("accumulator_precision must be Float32 or Float64") SatelliteGridding._spherical_ka_options(
            spec, backend, nothing, Float16)
        @test SatelliteGridding._spherical_ka_options(spec, backend, Float32, nothing).accumulator_precision === Float64
    end

    @testset "CUDA" begin
        available = try
            @eval begin
                using CUDA
                CUDA.functional()
            end
        catch
            false
        end
        if available
            @eval CUDA.allowscalar(false)
            backend = Base.invokelatest(resolve_backend, "cuda")
            Base.invokelatest(_test_healpix_ka_end_to_end, backend; compare_float32=true)
            for precision in (nothing, Float32)
                Base.invokelatest(_test_healpix_moment_precision, backend, precision)
                Base.invokelatest(_test_healpix_accumulator_cases, backend, precision)
                Base.invokelatest(_test_healpix_rolling_offset, backend, precision)
                Base.invokelatest(_test_healpix_rolling_mean_correction, backend, precision)
                Base.invokelatest(_test_healpix_tiny_weight, backend, precision)
            end
            Base.invokelatest(_healpix_test_precision, backend)
        else
            @test_skip available
        end
    end
end
