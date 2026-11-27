function _test_overlap_statistics(dir, backend; accumulator_precision=nothing)
    rect = GridSpec(lat_min=0.0, lat_max=2.0, lon_min=0.0, lon_max=2.0, dlat=1.0, dlon=1.0)
    cube = CubedSphereGridSpec(Nc=12, T=Float64)
    hp = HEALPixGridSpec(nside=8, T=Float64)
    hp_center = pix2ang_nest(8, 277)
    d1, d3, d4 = DateTime.(Date(2020, 1, 1) .+ Day.((0, 2, 3)))
    tiled = TimeSpec(d1, d3, Day(2))
    rolling = TimeSpec(d1, d4, Day(1); sample_step=Day(1), window_halfwidth_days=1)
    config = _healpix_ka_signal_config(dir, "stats")
    output = Dict{String,Any}("statistics" => ["sum_w2", "n_soundings"], "precision" => "float64")
    config.options["output"] = output
    path = joinpath(dir, "overlap.nc")
    for (spec, method, straddle, center) in (
        (rect, SubpixelGridding(), (0.8, 1.2, 0.8, 1.2), (0.5, 0.5)),
        (rect, CircularFootprintGridding(), (0.8, 1.2, 0.8, 1.2), (0.5, 0.5)),
        (cube, SubpixelGridding(), (-10.2, -9.8, -0.2, 0.2), (-10.1, 0.1)),
        (hp, SubpixelGridding(), (44.8, 45.2, -0.01, 0.25), hp_center),
    )
        lon, lat = center
        whole = (lon - 0.01, lon + 0.01, lat - 0.01, lat + 0.01)
        cell = if spec === rect
            1
        else
            S = SatelliteGridding
            S.locate(S._device_locator(spec), S._sphere_unit_vector(spec, lat, lon)...)
        end
        for (day, footprints, values) in ((1, [whole], [7.0]),
                                          (3, [straddle, straddle], [1.0, 3.0]),
                                          (4, [straddle, straddle], [5.0, 7.0]))
            _write_healpix_ka_signal_input(joinpath(dir, "stats_2020-01-0$day.nc"), footprints, values)
        end
        for (times, end_slice, compute_std) in ((tiled, 2, true), (rolling, 4, true), (rolling, 4, false))
            @testset "$(nameof(typeof(spec))) $(nameof(typeof(method))) std=$compute_std slices=$end_slice" begin
                grid_l2(config, spec, times; outfile=path, n_oversample=2, footprint_method=method,
                        backend, accumulator_precision, compute_std, progress=false, verbose=false)
                Dataset(path) do ds
                    read_slice(name, t) = vec(ds[name][ntuple(_ -> Colon(), ndims(ds[name]) - 1)..., t])
                    @test eltype(ds["sum_w2"].var) === Float64
                    @test eltype(ds["n_soundings"].var) === Int32
                    @test occursin("Kish", ds["sum_w2"].attrib["comment"])
                    @test occursin("sub-samples", ds["n_soundings"].attrib["comment"])
                    @test ds.attrib["accumulator_precision"] == lowercase(string(something(accumulator_precision, Float64)))
                    if spec === rect
                        @test !haskey(ds.attrib, "geometry_precision")
                    else
                        @test ds.attrib["geometry_precision"] == "float64"
                    end
                    @test sum(read_slice("n", 1)) == sum(read_slice("sum_w2", 1)) ==
                          sum(read_slice("n_soundings", 1)) == 1
                    @test read_slice("n", 1)[cell] == read_slice("sum_w2", 1)[cell] ==
                          read_slice("n_soundings", 1)[cell] == 1
                    weights = read_slice("n", end_slice)
                    squares = read_slice("sum_w2", end_slice)
                    counts = read_slice("n_soundings", end_slice)
                    @test weights[cell] == 1
                    @test squares[cell] == 1/4
                    @test weights[cell]^2 / squares[cell] == 4
                    @test counts[cell] == 4
                    @test all(iszero, squares[iszero.(weights)])
                    @test all(iszero, counts[iszero.(weights)])
                    @test read_slice("signal", end_slice)[cell] ≈ 4.0
                    if compute_std
                        @test read_slice("signal_std", end_slice)[cell] ≈ sqrt(5.0) rtol=1e-6
                    end
                end
            end
        end
    end
end

@kernel function _increment_test_soundings!(counts)
    SatelliteGridding._increment_soundings!(counts, 1, Val(true))
end

@testset "Sounding count overflow" begin
    S = SatelliteGridding
    max_count = typemax(Int32)
    overflow = OverflowError("n_soundings exceeds the Int32 limit ($max_count)")
    counts = [max_count - Int32(1)]
    S._increment_soundings!(counts, 1, Val(false))
    @test counts == [max_count]
    @test S._check_overlap_slice((n_soundings=counts,)) === nothing
    S._increment_soundings!(counts, 1, Val(false))
    @test_throws overflow S._check_overlap_slice((n_soundings=counts,))

    backend = CPU()
    for start in (max_count - Int32(1), Int32(-1))
        counts[1] = start
        _increment_test_soundings!(backend)(counts; ndrange=4)
        KernelAbstractions.synchronize(backend)
        @test_throws overflow S._check_overlap_slice((n_soundings=counts,))
    end

    counts[1] = max_count
    window = (n_soundings=counts,)
    S._merge_overlap!(window, (n_soundings=Int32[1],))
    @test_throws overflow S._check_overlap_slice(window)
end

@testset "Opt-in footprint statistics" begin
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        for (backend, precision) in ((nothing, nothing), (CPU(), nothing), (CPU(), Float32))
            @testset "$backend $precision" begin
                _test_overlap_statistics(dir, backend; accumulator_precision=precision)
            end
        end
        @test_throws ArgumentError NetCDFOptions(; statistics=["n_soundings", "unknown"])
        @test_throws ArgumentError NetCDFOptions(; statistics="sum_w2")
        S = SatelliteGridding
        config = _healpix_ka_signal_config(dir, "stats")
        day = DateTime(2020, 1, 1)
        spec = HEALPixGridSpec(nside=8)
        times = TimeSpec(day, day, Day(1))
        for (name, T) in (("sum_w2", Float32), ("n_soundings", Int32))
            config.options["output"] = Dict("statistics" => [name])
            path = joinpath(dir, "single.nc")
            grid_l2(config, spec, times; outfile=path, backend=CPU(), progress=false, verbose=false)
            Dataset(path) do ds
                @test [key for key in keys(ds) if key in ("sum_w2", "n_soundings")] == [name]
                @test eltype(ds[name].var) === T
            end
            @test_throws ArgumentError S._validate_zarr_output(GridSpec(), S._output_options(config.options))
        end
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
            Base.invokelatest(_test_overlap_statistics, dir, backend)
        else
            @test_skip available
        end
    end
end
