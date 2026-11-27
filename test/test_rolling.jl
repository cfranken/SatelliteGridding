@testset "Rolling mean — centered window with day cache" begin

    @testset "TimeSpec rolling construction" begin
        tiled = TimeSpec(DateTime("2020-01-01"), DateTime("2020-03-01"), Dates.Day(8))
        @test !SatelliteGridding.is_rolling(tiled)

        roll = TimeSpec(DateTime("2020-01-01"), DateTime("2020-03-01"), Dates.Day(8);
                        sample_step=Dates.Day(1), window_halfwidth_days=7)
        @test SatelliteGridding.is_rolling(roll)
        @test roll.sample_step == Dates.Day(1)
        @test roll.window_halfwidth_days == 7

        # Exactly one of the two rolling fields → error
        @test_throws ErrorException TimeSpec(DateTime("2020-01-01"), DateTime("2020-03-01"),
                                             Dates.Day(8); sample_step=Dates.Day(1))
        # Negative half-width → error
        @test_throws ErrorException TimeSpec(DateTime("2020-01-01"), DateTime("2020-03-01"),
                                             Dates.Day(8); sample_step=Dates.Day(1),
                                             window_halfwidth_days=-1)
    end

    @testset "moment helpers vs brute force" begin
        # One cell / one variable, weighted samples spread over 3 days.
        days = [(x=Float32[1, 2, 3], w=Float32[1, 1, 2]),
                (x=Float32[5, 4],    w=Float32[3, 1]),
                (x=Float32[0, 10, 2], w=Float32[1, 1, 1])]

        S = SatelliteGridding
        spec = GridSpec(lat_min=0f0, lat_max=1f0, lon_min=0f0, lon_max=1f0, dlat=1f0, dlon=1f0)
        acc = S.new_accumulator(spec, nothing, nothing, 1, true)
        window = S._window_moments(acc)
        S.reset!(window)
        @test eltype(window.weights) == eltype(window.mean) == eltype(window.m2) == Float64
        for day in days
            S.reset!(acc)
            vals = reshape(day.x, :, 1)
            for i in eachindex(day.x)
                S._welford_update!(acc.data, acc.m2, acc.weights, true, 1, 1, vals, i, 1, day.w[i])
            end
            S._merge_day!(window, S.day_moments(acc))
        end
        means, stds, weights = S.finish!(window, acc)

        allx = vcat((d.x for d in days)...); allw = vcat((d.w for d in days)...)
        W = sum(allw); mu = sum(allw .* allx) / W
        sd = sqrt(sum(allw .* (allx .- mu) .^ 2) / W)

        @test weights[1, 1] == W
        @test means[1, 1, 1] ≈ mu atol = 1e-4
        @test stds[1, 1, 1] ≈ sd atol = 1e-4
    end

    @testset "End-to-end rolling pipeline" begin
        test_dir = mktempdir()
        # One footprint per day in the same 2°×2° cell; signal = 10·day.
        for d in 1:5
            f = joinpath(test_dir, "L2_2020-01-0$(d)_x.nc")
            ds = Dataset(f, "c")
            defDim(ds, "x", 1); defDim(ds, "nv", 4)
            vlat = defVar(ds, "lat", Float32, ("x",)); vlon = defVar(ds, "lon", Float32, ("x",))
            vlb = defVar(ds, "lat_bnds", Float32, ("x", "nv")); vob = defVar(ds, "lon_bnds", Float32, ("x", "nv"))
            vs = defVar(ds, "sif", Float32, ("x",))
            vlat[:] = Float32[1.0]; vlon[:] = Float32[1.0]
            vlb[:, :] = Float32[0.75 0.75 1.25 1.25]; vob[:, :] = Float32[0.75 1.25 1.25 0.75]
            vs[:] = Float32[10.0 * d]
            close(ds)
        end

        cfg = joinpath(test_dir, "c.toml")
        write(cfg, """
        filePattern = "L2_YYYY-MM-DD_*.nc"
        folder = "$(replace(test_dir, "\\" => "/"))/"
        [basic]
        lat = "lat"
        lon = "lon"
        lat_bnd = "lat_bnds"
        lon_bnd = "lon_bnds"
        [grid]
        sif = "sif"
        """)

        config = load_config(cfg)
        gs = GridSpec(lat_min=0f0, lat_max=2f0, lon_min=0f0, lon_max=2f0, dlat=2f0, dlon=2f0)
        ts = TimeSpec(DateTime("2020-01-01"), DateTime("2020-01-05"), Dates.Day(1);
                      sample_step=Dates.Day(1), window_halfwidth_days=1)

        out = joinpath(test_dir, "roll.nc")
        grid(config, gs, ts, SubpixelGridding(n_oversample=4); compute_std=true, outfile=out)

        ds = Dataset(out)
        sif = ds["sif"][1, 1, :]
        n = ds["n"][1, 1, :]
        t = ds["time"][:]; ws = ds["window_start"][:]; we = ds["window_end"][:]

        # Centered ±1-day windows clip at the data edges (days 0 and 6 are absent).
        @test length(t) == 5
        @test sif ≈ Float32[15, 20, 30, 40, 45] atol = 1e-3
        @test n ≈ Float32[2, 3, 3, 3, 2]

        # time = window center; window_start/window_end bracket it symmetrically.
        @test t[3] == DateTime("2020-01-03")
        @test ws[3] == DateTime("2020-01-02")
        @test we[3] == DateTime("2020-01-04")
        @test ws[1] == DateTime("2019-12-31")   # extends before start_date
        @test we[5] == DateTime("2020-01-06")   # extends past stop_date

        # std at day 3 over {20,30,40}: population std = sqrt(200/3).
        @test ds["sif_std"][1, 1, 3] ≈ sqrt(200 / 3) atol = 1e-2
        close(ds)

        # KA mean-only rolling matches the sequential means.
        outka = joinpath(test_dir, "roll_ka.nc")
        grid(config, gs, ts, SubpixelGridding(n_oversample=4); backend=CPU(), outfile=outka)
        dk = Dataset(outka)
        @test dk["sif"][1, 1, :] ≈ sif atol = 1e-3
        close(dk)

        grid(config, gs, ts, SubpixelGridding(n_oversample=4);
             backend=CPU(), compute_std=true, outfile=outka)
        Dataset(outka) do dk
            @test dk["sif"][1, 1, :] ≈ sif atol = 1e-3
            @test dk["n"][1, 1, :] == n
            @test dk["sif_std"][1, 1, 3] ≈ sqrt(200 / 3) atol = 1e-4
        end

        rm(test_dir; recursive=true)
    end
end

function _weighted_population_reference(values, weights)
    # Evaluate against the given observation values, including their input rounding.
    x, w = BigFloat.(values), BigFloat.(weights)
    total = sum(w)
    mean = sum(w .* x) / total
    std = sqrt(sum(w .* (x .- mean) .^ 2) / total)
    Float64(mean), Float64(std)
end

function _test_rectangular_l2_std(dir, backend, accumulator_precision)
    spec = GridSpec(lat_min=0f0, lat_max=7f0, lon_min=0f0, lon_max=7f0, dlat=1f0, dlon=1f0)
    config = _healpix_ka_signal_config(dir, "rect_std")
    target = (0.25, 0.75, 0.25, 0.75)
    half_target = (0.65, 1.35, 0.25, 0.75)  # Half in each of two cells.
    singleton = (3.65, 4.35, 0.25, 0.75)   # One sounding, two fractional cells.
    wide = (2.1, 5.9, 2.1, 5.9)           # Sixteen cells, one sub-sample per cell.
    dates = DateTime(2020, 1, 1):Day(1):DateTime(2020, 1, 4)
    tiled = TimeSpec(first(dates), first(dates), Day(4))
    rolling = TimeSpec(dates[2], dates[3], Day(1); sample_step=Day(1), window_halfwidth_days=1)
    for (offset, spread) in ((410.0, 0.05), (410.0, 3.0), (1e8, 0.05))
        offset == 1e8 && accumulator_precision === Float32 && continue
        # Different daily means, unequal numbers of soundings and overlap weights.
        # Test Float32 accumulators against observations stored in Float32.
        T = something(accumulator_precision, Float64)
        daily = [T.(offset .+ spread .* range(-0.9 + day/5, 0.6 + day/5; length=day+3))
                 for day in 1:4]
        for (day, date) in enumerate(dates)
            values = daily[day]
            rectangles = [fill(target, length(values) - 1); [half_target]]
            if day == 2
                append!(rectangles, (singleton, wide))
                values = [values; T(offset); T(offset)]
            end
            _write_healpix_ka_signal_input(joinpath(dir, "rect_std_$(Date(date)).nc"), rectangles, values)
        end
        config.options["output"] = Dict("precision" => "float64",
            "significant_digits_mean" => "none", "significant_digits_std" => "none")
        for (times, windows) in ((tiled, (1:4,)), (rolling, (1:3, 2:4)))
            ka, seq = joinpath(dir, "std_ka.nc"), joinpath(dir, "std_seq.nc")
            options = (; compute_std=true, n_oversample=4, progress=false, verbose=false)
            grid_l2(config, spec, times; backend, accumulator_precision, outfile=ka, options...)
            # Tiled Float32 Welford cannot resolve a spread of 0.05 at 1e8.
            grid_l2(config, spec, times; outfile=seq, options...)
            Dataset(ka) do dk
                Dataset(seq) do ds
                    @test dk["n"][:, :, :] == ds["n"][:, :, :]
                    for (time, days) in enumerate(windows)
                        values = vcat(daily[days]...)
                        weights = vcat(([ones(length(day) - 1); 0.5] for day in daily[days])...)
                        mean, std = _weighted_population_reference(values, weights)
                        @test dk["n"][1, 1, time] == sum(weights)
                        @test dk["signal"][1, 1, time] ≈ mean rtol=2e-7
                        @test dk["signal_std"][1, 1, time] ≈ std rtol=1e-4
                        if times === rolling
                            @test ds["signal_std"][1, 1, time] ≈ std rtol=1e-4
                        end
                        if times === rolling || offset == 410.0
                            populated = dk["n"][:, :, time] .> 0
                            @test all(isapprox.(dk["signal"][:, :, time][populated],
                                                ds["signal"][:, :, time][populated]; rtol=4e-7))
                            @test all(isapprox.(dk["signal_std"][:, :, time][populated],
                                                ds["signal_std"][:, :, time][populated]; rtol=3e-3, atol=1e-7))
                        end
                        for col in (4, 5)
                            @test dk["n"][col, 1, time] == 0.5
                            @test dk["signal_std"][col, 1, time] == ds["signal_std"][col, 1, time] == 0
                        end
                        @test dk["n"][3:6, 3:6, time] == fill(1/16, 4, 4)
                        @test all(iszero, dk["signal_std"][3:6, 3:6, time])
                    end
                end
            end
        end
    end
end

function _test_l2_circular_weight(dir, backend)
    spec = GridSpec(lat_min=0f0, lat_max=2f0, lon_min=0f0, lon_max=2f0, dlat=1f0, dlon=1f0)
    config = _healpix_ka_signal_config(dir, "circle")
    rect = (0.3, 0.7, 0.3, 0.7)
    _write_healpix_ka_signal_input(joinpath(dir, "circle_2020-01-01.nc"), fill(rect, 1000), fill(410.1f0, 1000))
    date = DateTime(2020, 1, 1)
    times = TimeSpec(date, date, Day(1))
    for b in (nothing, backend), compute_std in (false, true)
        output = joinpath(dir, "circle_out.nc")
        grid_l2(config, spec, times; backend=b, compute_std, outfile=output,
                footprint_method=CircularFootprintGridding(20), verbose=false, progress=false)
        Dataset(output) do ds
            @test ds["n"][:, :, 1] == Float32[1000 0; 0 0]
            @test ds["signal"][1, 1, 1] == 410.1f0
            compute_std && @test ds["signal_std"][1, 1, 1] == 0
        end
    end
end

function _test_ka_mean_refinement(dir, backend)
    # The rolling representation must preserve a refined mean at the exact
    # lower bound, even when rebasing would round away the reference.
    S = SatelliteGridding
    mean = fill(3f0, 1, 1, 1)
    residual = fill(-2f0, size(mean))
    weights = ones(Float32, 1, 1)
    reference = Float32[-eps(1f0) / 8]
    bounds = (lo=fill(S._ka_ordered_key(1f0), size(mean)),
              hi=fill(S._ka_ordered_key(3f0), size(mean)))
    S._correct_ka_moments_ka!(backend, 1)(mean, nothing, residual, weights, reference,
                                       bounds, Val(true); ndrange=1)
    KernelAbstractions.synchronize(backend)
    @test Float64(mean[1]) + Float64(residual[1]) == 1.0

    start = DateTime(2020, 1, 1)
    next = start + Day(1)
    tiled = TimeSpec(start, start, Day(2))
    rolling = TimeSpec(start, next, Day(1); sample_step=Day(1), window_halfwidth_days=0)
    rect = GridSpec(lat_min=0f0, lat_max=1f0, lon_min=0f0, lon_max=2f0, dlat=1f0, dlon=1f0)
    healpix = HEALPixGridSpec(nside=8)
    for (spec, anchor, target, cell) in (
        (rect, (0.2, 0.7, 0.2, 0.7), (1.2, 1.7, 0.2, 0.7), (2, 1)),
        (healpix, _healpix_cell_rectangle(8, 500), _healpix_cell_rectangle(8, 100), (101,)),
    )
        config = _healpix_ka_signal_config(dir, "refinement")
        config.options["output"] = Dict("precision" => "float64",
            "significant_digits_mean" => "none", "significant_digits_std" => "none")
        _write_healpix_ka_signal_input(joinpath(dir, "refinement_$(Date(start)).nc"), [anchor], [410.1f0])
        target_path = joinpath(dir, "refinement_$(Date(next)).nc")
        targets = fill(target, 100_000)
        values = fill(1f0, 100_000)
        _write_healpix_ka_signal_input(target_path, targets, values)
        outfile = joinpath(dir, "refinement_out.nc")
        for (times, time_index) in ((tiled, 1), (rolling, 2))
            grid_l2(config, spec, times; backend, accumulator_precision=Float32, compute_std=false,
                    n_oversample=1, outfile, verbose=false, progress=false)
            Dataset(outfile) do ds
                @test ds["n"][cell..., time_index] == 100_000
                @test abs(ds["signal"][cell..., time_index] - 1f0) <= eps(1f0)
                @test !haskey(ds, "signal_std")
            end
        end

        values[end] = 0.5f0
        rm(target_path)
        _write_healpix_ka_signal_input(target_path, targets, values)
        for (times, time_index) in ((tiled, 1), (rolling, 2)),
                (accumulator_precision, compute_std) in ((Float32, false), (Float32, true), (Float64, true))
            grid_l2(config, spec, times; backend, accumulator_precision, compute_std,
                    n_oversample=1, outfile, verbose=false, progress=false)
            Dataset(outfile) do ds
                @test ds["n"][cell..., time_index] == 100_000
                observed_mean = ds["signal"][cell..., time_index]
                if accumulator_precision === Float32
                    @test 0.5 <= observed_mean <= 1
                else
                    @test observed_mean ≈ 0.999995 atol=1e-12
                    @test ds["signal_std"][cell..., time_index] ≈ 0.001581130924 atol=1e-12
                end
            end
        end
    end
end

function _test_rolling_selection_without_std(dir)
    date = DateTime(2020, 1, 1)
    path = joinpath(dir, "selection_$(Date(date)).nc")
    _write_healpix_ka_signal_input(path, fill((0.2, 0.7, 0.2, 0.7), 2), [10.0, 30.0])
    Dataset(path, "a") do ds
        defVar(ds, "wide", Float64, ("sounding",))[:] = [1.0, 1e40]
    end
    config = _healpix_ka_signal_config(dir, "selection")
    config.grid_vars["wide"] = "wide"
    config.options["output"] = Dict("precision" => "float64", "significant_digits_mean" => "none")
    spec = GridSpec(lat_min=0f0, lat_max=1f0, lon_min=0f0, lon_max=1f0, dlat=1f0, dlon=1f0)
    times = TimeSpec(date, date, Day(1); sample_step=Day(1), window_halfwidth_days=0)
    outfile = joinpath(dir, "selection_out.nc")
    for compute_std in (false, true)
        grid_l2(config, spec, times; backend=nothing, compute_std, n_oversample=1,
                outfile, verbose=false, progress=false)
        Dataset(outfile) do ds
            @test ds["n"][1, 1, 1] == 2
            @test ds["signal"][1, 1, 1] == 20
            @test ds["wide"][1, 1, 1] == 5e39
        end
    end
end

function _test_spherical_rolling_std(dir, backend)
    spec = HEALPixGridSpec(nside=8)
    config = _healpix_ka_signal_config(dir, "spherical_std")
    target = _healpix_cell_rectangle(8, 282)
    singleton = _healpix_cell_rectangle(8, 100)
    # At n=4, this grazing footprint has six of 16 samples in NESTED ID 282.
    grazing = (-0.03, 0.04, -0.025, 0.02)
    days = DateTime(2020, 1, 1):Day(1):DateTime(2020, 1, 3)
    times = TimeSpec(days[2], days[2], Day(1); sample_step=Day(1), window_halfwidth_days=1)
    config.options["output"] = Dict("precision" => "float64",
        "significant_digits_mean" => "none", "significant_digits_std" => "none")
    for (offset, spread) in ((410.0, 0.05), (410.0, 3.0), (1e8, 0.05))
        deviations = [-0.7 0.1 0.5; -0.8 0.3 0.9]
        values = offset .+ spread .* deviations
        for (idx, date) in enumerate(days)
            rectangles = idx == 2 ? [target, grazing, singleton] : [target, grazing]
            signal = idx == 2 ? [values[:, idx]; offset] : values[:, idx]
            _write_healpix_ka_signal_input(joinpath(dir, "spherical_std_$(Date(date)).nc"), rectangles, signal)
        end
        weights = repeat([1.0, 6 / 16], length(days))
        expected_mean, expected_std = _weighted_population_reference(vec(values), weights)
        output = joinpath(dir, "spherical_out.nc")
        grid_l2(config, spec, times; backend, compute_std=true, n_oversample=4, outfile=output,
                verbose=false, progress=false)
        Dataset(output) do ds
            @test ds["n"][283, 1] == sum(weights)
            @test ds["signal"][283, 1] ≈ expected_mean rtol=1e-12
            @test ds["signal_std"][283, 1] ≈ expected_std rtol=1e-4
            @test ds["n"][101, 1] == 1
            @test ds["signal_std"][101, 1] == 0
        end
    end
end

@testset "KA Level-2 centered statistics" begin
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        backend = resolve_backend("cpu")
        for precision in (nothing, Float32)
            _test_rectangular_l2_std(dir, backend, precision)
        end
        _test_l2_circular_weight(dir, backend)
        _test_ka_mean_refinement(dir, backend)
        _test_rolling_selection_without_std(dir)
        _test_spherical_rolling_std(dir, backend)

        if isdefined(Main, :CUDA) && Main.CUDA.functional()
            cuda = Base.invokelatest(resolve_backend, "cuda")
            Base.invokelatest(_test_rectangular_l2_std, dir, cuda, nothing)
            Base.invokelatest(_test_l2_circular_weight, dir, cuda)
            Base.invokelatest(_test_spherical_rolling_std, dir, cuda)
        end
    end
end
