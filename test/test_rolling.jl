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

        function day_moments(d)
            mean = zeros(Float32, 1, 1, 1); M2 = zeros(Float32, 1, 1, 1); w = zeros(Float32, 1, 1)
            vals = reshape(d.x, :, 1)
            for i in eachindex(d.x)
                SatelliteGridding._welford_update!(mean, M2, w, true, 1, 1, vals, i, 1, d.w[i])
            end
            S1 = Array{Float32}(undef, 1, 1, 1); S2 = Array{Float32}(undef, 1, 1, 1)
            SatelliteGridding.welford_to_moments!(S1, S2, mean, M2, w, true)
            (S0=copy(w), S1=S1, S2=S2)
        end

        ms = day_moments.(days)
        S0 = sum(m.S0 for m in ms); S1 = sum(m.S1 for m in ms); S2 = sum(m.S2 for m in ms)
        SatelliteGridding.finalize_moments!(S1, S2, S0, true)

        allx = vcat((d.x for d in days)...); allw = vcat((d.w for d in days)...)
        W = sum(allw); mu = sum(allw .* allx) / W
        sd = sqrt(sum(allw .* (allx .- mu) .^ 2) / W)

        @test S1[1, 1, 1] ≈ mu atol = 1e-4
        @test S2[1, 1, 1] ≈ sd atol = 1e-4
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

        # KA + std + rolling is unsupported and must error clearly.
        @test_throws ErrorException grid(config, gs, ts, SubpixelGridding(n_oversample=4);
                                         backend=CPU(), compute_std=true,
                                         outfile=joinpath(test_dir, "err.nc"))

        rm(test_dir; recursive=true)
    end
end
