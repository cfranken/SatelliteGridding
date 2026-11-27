@testset "Coarsen a NESTED HEALPix archive" begin
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        _write_l2_hash_fixture(dir)
        vars = OrderedDict("signal" => "signal", "offset" => "offset")
        cfg = DataSourceConfig(Dict("lat_bnd" => "lat", "lon_bnd" => "lon"), vars,
                               [FilterRule("qa", :gt, 0.8)], "input_YYYY-MM-DD.nc", dir)
        output = Dict{String,Any}("precision" => "float64", "significant_digits_mean" => "none",
            "significant_digits_std" => "none", "fill_value" => NaN,
            "statistics" => ["sum_w2", "n_soundings"])
        cfg.options["output"] = output
        times = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 5), Day(1))
        opts = (; compute_std=true, n_oversample=7, progress=false, verbose=false)
        fine, coarse, direct = joinpath.(dir, ("fine.nc", "coarse.nc", "direct.nc"))
        for backend in (nothing, CPU())
            grid_l2(cfg, HEALPixGridSpec(nside=16, T=Float64, cell_bounds=true), times;
                    outfile=fine, backend, opts...)
            for (levels, nside) in ((1, 8), (2, 4))
                if levels == 1
                    @test_logs coarsen_healpix(fine, coarse)
                else
                    @test_logs ProgressTestCLI.main(
                        ["coarsen", "--input", fine, "--output", coarse, "--nside", string(nside)])
                end
                grid_l2(cfg, HEALPixGridSpec(; nside, T=Float64, cell_bounds=true), times;
                        outfile=direct, backend, opts...)
                Dataset(coarse) do actual
                    Dataset(direct) do expected
                        @test actual.attrib["coarsen_source"] == abspath(fine)
                        @test actual.attrib["coarsen_levels"] == levels
                        @test actual.attrib["coarsen_source_nside"] == 16
                        @test actual.attrib["coarsen_target_nside"] == nside
                        @test actual.attrib["source_versions"] == expected.attrib["source_versions"]
                        @test haskey(actual.attrib, "coarsen_command") == (levels == 2)
                        for name in ("nside", "healpix_level", "healpix_ordering", "Conventions")
                            @test actual.attrib[name] == expected.attrib[name]
                        end
                        for name in ("time", "window_start", "window_end", "cell_ids",
                                     "lon", "lat", "lon_bounds", "lat_bounds")
                            @test actual[name].var[:] == expected[name].var[:]
                            @test Dict(actual[name].attrib) == Dict(expected[name].attrib)
                        end
                        @test actual["cell_area"][] == expected["cell_area"][]
                        @test Dict(actual["healpix"].attrib) == Dict(expected["healpix"].attrib)
                        weights = expected["n"][:, :]
                        empty = iszero.(weights)
                        @test actual["n"][:, :] ≈ weights rtol=1e-15
                        @test all(iszero, actual["n"][:, 4:5])
                        @test all(iszero, actual["sum_w2"][:, :][empty])
                        @test all(iszero, actual["n_soundings"][:, :][empty])
                        @test all(actual["sum_w2"][:, :][.!empty] .<= expected["sum_w2"][:, :][.!empty] .+ 1e-15)
                        @test all(actual["n_soundings"][:, :] .>= expected["n_soundings"][:, :])
                        @test occursin("Underestimates", actual["sum_w2"].attrib["comment"])
                        @test occursin("Overcounts", actual["n_soundings"].attrib["comment"])
                        for (name, role) in (("n", "weight_sum"), ("sum_w2", "sum_w2"),
                                             ("n_soundings", "n_soundings"), ("signal_std", "population_std"))
                            @test actual[name].attrib["satellitegridding_role"] == role
                            @test expected[name].attrib["satellitegridding_role"] == role
                        end
                        for name in ("sum_w2", "n_soundings")
                            @test occursin("standard error", actual[name].attrib["coarsening_approximation"])
                            @test !haskey(expected[name].attrib, "coarsening_approximation")
                        end
                        @test actual["signal_std"].attrib["satellitegridding_mean"] == "signal"
                        for name in ("signal", "offset", "signal_std", "offset_std")
                            @test eltype(actual[name].var) === Float64
                            @test all(ismissing, actual[name][:, :][empty])
                            @test actual[name].attrib["coordinates"] == expected[name].attrib["coordinates"]
                            @test get(actual[name].attrib, "units", nothing) == get(expected[name].attrib, "units", nothing)
                            @test all(isapprox.(actual[name][:, :][.!empty], expected[name][:, :][.!empty];
                                                rtol=1e-12, atol=name == "offset_std" ? 4e-12 : 1e-15))
                        end
                    end
                end
            end
        end
        @test_throws ArgumentError coarsen_healpix(fine, fine)
        @test_throws ArgumentError coarsen_healpix(fine, coarse; levels=2, nside=4)
        @test_throws ArgumentError coarsen_healpix(fine, coarse; nside=16)
        @test_throws ArgumentError coarsen_healpix(fine, coarse; levels=5)

        # Filled science on a positive-weight child preserves n and a valid
        # mean when only std is missing.
        parent = Dataset(fine, "a") do ds
            child = findfirst(>(0), ds["n"][:, 1])
            ds["signal_std"].var[child, 1] = NaN
            ds["offset"].var[child, 1] = NaN
            for name in ("sum_w2", "n_soundings")
                ds[name].var[child, 1] = ds[name].attrib["_FillValue"]
            end
            cld(child, 16)
        end
        coarsen_healpix(fine, coarse; levels=2)
        Dataset(coarse) do ds
            @test ds["n"][parent, 1] > 0
            @test isfinite(ds["signal"][parent, 1])
            @test ismissing(ds["signal_std"][parent, 1])
            @test ismissing(ds["offset"][parent, 1])
            @test ismissing(ds["offset_std"][parent, 1])
            @test ismissing(ds["sum_w2"][parent, 1])
            @test ismissing(ds["n_soundings"][parent, 1])
        end

        output["precision"] = "float32"
        output["significant_digits_mean"] = 4
        output["significant_digits_std"] = 3
        grid_l2(cfg, HEALPixGridSpec(nside=16), times; outfile=fine, opts...)
        warning = r"Float32 output, Float32 accumulation, decimal rounding.*exact-archive"
        @test_logs (:warn, warning) coarsen_healpix(fine, coarse; levels=1)
        @test_logs (:warn, warning) coarsen_healpix(coarse, direct; levels=1)
        Dataset(direct) do ds
            @test ds.attrib["coarsen_source"] == abspath(fine)
            @test ds.attrib["coarsen_levels"] == 2
            @test ds.attrib["coarsen_source_nside"] == 16
            @test ds.attrib["coarsen_target_nside"] == ds.attrib["nside"] == 4
        end
    end
end

@testset "Coarsening precision and science fields named like statistics" begin
    mktempdir(joinpath(@__DIR__, "..", ".scratch")) do dir
        path = joinpath(dir, "precise_2020-01-01.nc")
        footprints = map(280:282) do cell
            lon, lat = pix2ang_nest(8, cell)
            (lon - 0.01, lon + 0.01, lat - 0.01, lat + 0.01)
        end
        _write_healpix_ka_signal_input(path, footprints, Float64(2^52) .+ (0:2))
        Dataset(path, "a") do ds
            defVar(ds, "marker", Float64, ("sounding",))[:] = [-999.0, -998.0, -997.0]
            ds["marker"].attrib["long_name"] = "Number of footprints touching the grid cell"
        end
        cfg = _healpix_ka_signal_config(dir, "precise")
        cfg.grid_vars["n_soundings"] = "marker"
        cfg.options["metadata"] = Dict("variable_attributes" => ["long_name", "comment"])
        cfg.options["output"] = Dict("precision" => "float64", "significant_digits_mean" => "none",
                                     "significant_digits_std" => "none", "fill_value" => NaN)
        day = DateTime(2020, 1, 1)
        times = TimeSpec(day, day, Day(1))
        fine, coarse, direct = joinpath.(dir, ("fine.nc", "coarse.nc", "direct.nc"))
        opts = (; compute_std=true, backend=CPU(), n_oversample=4, progress=false, verbose=false)
        grid_l2(cfg, HEALPixGridSpec(nside=8), times; outfile=fine, opts...)
        grid_l2(cfg, HEALPixGridSpec(nside=4), times; outfile=direct, opts...)
        @test_logs coarsen_healpix(fine, coarse)
        parent = 71
        Dataset(fine) do ds
            @test ds["signal"][281:283, 1] == Float64(2^52) .+ (0:2)
            @test ds["n_soundings"][281, 1] == -999.0
            @test ds["n_soundings"].attrib["satellitegridding_role"] == "mean"
            @test ds["n_soundings"].attrib["satellitegridding_std"] == "n_soundings_std"
        end
        Dataset(coarse) do actual
            Dataset(direct) do expected
                for name in ("signal", "n_soundings", "n_soundings_std", "n")
                    @test actual[name][parent, 1] == expected[name][parent, 1]
                end
                @test actual["signal_std"][parent, 1] ≈ expected["signal_std"][parent, 1] ≈ sqrt(2/3)
                @test eltype(actual["n_soundings"].var) === Float64
                @test actual["n_soundings"][parent, 1] == -998.0
            end
        end
        grid_l2(cfg, HEALPixGridSpec(nside=8), times; outfile=fine, opts..., compute_std=false)
        @test_logs coarsen_healpix(fine, coarse)
        Dataset(coarse) do ds
            @test ds["n_soundings"][parent, 1] == -998.0
            @test ds["n_soundings"].attrib["satellitegridding_role"] == "mean"
            @test !haskey(ds["n_soundings"].attrib, "satellitegridding_std")
        end

        config = joinpath(dir, "precise.toml")
        write(config, """
            filePattern = "precise_YYYY-MM-DD.nc"
            folder = "$dir"
            [basic]
            lat_bnd = "lat_bounds"
            lon_bnd = "lon_bounds"
            [grid]
            signal = "signal"
            [output]
            precision = "float64"
            significant_digits_mean = "none"
            significant_digits_std = "none"
            fill_value = nan
            """)
        args = ["l2", "--config", config, "--outFile", fine, "--gridType", "healpix", "--nside", "8",
                "--startDate", "2020-01-01", "--stopDate", "2020-01-01", "--nOversample", "4", "--quiet"]
        ProgressTestCLI.main(args)
        Dataset(fine) do ds
            @test ds.attrib["accumulator_precision"] == "float32"
            @test ds.attrib["geometry_precision"] == "float64"
            @test eltype(ds["n"].var) === Float64
        end
        @test_logs (:warn, r"Float32 accumulation\. See the exact-archive") coarsen_healpix(fine, coarse)
        ProgressTestCLI.main([args; "--accumulator-precision"; "float64"])
        Dataset(fine) do ds
            @test ds.attrib["accumulator_precision"] == ds.attrib["geometry_precision"] == "float64"
            @test ds["signal"][281:283, 1] == Float64(2^52) .+ (0:2)
        end
        @test_logs coarsen_healpix(fine, coarse)
        rolling = [args; "--sample_step"; "1"]
        @test_throws ArgumentError("--accumulator-precision float32 is unsupported for sequential rolling, which accumulates daily moments in Float64") ProgressTestCLI.main(
            [rolling; "--backend"; "sequential"; "--accumulator-precision"; "float32"])
        ProgressTestCLI.main([rolling; "--accumulator-precision"; "float64"])
        Dataset(fine) do ds
            @test ds.attrib["accumulator_precision"] == "float64"
            @test ds["signal"][281:283, 1] == Float64(2^52) .+ (0:2)
        end
    end
end

@testset "Explicit and legacy coarsening pairs" begin
    mktempdir(joinpath(@__DIR__, "..", ".scratch")) do dir
        fine, reordered, coarse, invalid = joinpath.(dir, ("fine.nc", "reordered.nc", "coarse.nc", "invalid.nc"))
        direct = joinpath(dir, "direct.nc")
        path = joinpath(dir, "roles_2020-01-01.nc")
        footprints = map(280:281) do cell
            lon, lat = pix2ang_nest(8, cell)
            (lon - 0.01, lon + 0.01, lat - 0.01, lat + 0.01)
        end
        _write_healpix_ka_signal_input(path, footprints, [2.0, 4.0])
        copied_comment = "sqrt(M2/n), where M2 is the weighted sum of squared deviations and n is the sum of weights. "
        Dataset(path, "a") do ds
            defVar(ds, "error", Float64, ("sounding",); attrib=["comment" => copied_comment])[:] = [10.0, 10.0]
        end
        cfg = _healpix_ka_signal_config(dir, "roles")
        cfg.grid_vars["signal_error"] = "error"
        cfg.options["metadata"] = Dict("variable_attributes" => ["comment"])
        cfg.options["output"] = Dict("precision" => "float64", "significant_digits_mean" => "none",
                                     "significant_digits_std" => "none", "fill_value" => NaN)
        day = DateTime(2020, 1, 1)
        spec, times = HEALPixGridSpec(nside=8, T=Float64), TimeSpec(day, day, Day(1))
        opts = (; backend=CPU(), n_oversample=4, progress=false, verbose=false)
        parent_spec = HEALPixGridSpec(nside=4, T=Float64)
        grid_l2(cfg, spec, times; outfile=fine, compute_std=true, opts...)
        grid_l2(cfg, parent_spec, times; outfile=direct, compute_std=true, opts...)
        Dataset(fine) do src
            Dataset(reordered, "c"; attrib=Dict(src.attrib)) do dst
                for name in sort(collect(keys(src)))
                    defVar(dst, src[name])
                end
            end
        end
        for tagged in (true, false)
            @test_logs coarsen_healpix(reordered, coarse; legacy_std=tagged ? :none : :by_name)
            Dataset(coarse) do ds
                Dataset(direct) do expected
                    for name in ("n", "signal", "signal_std", "signal_error", "signal_error_std")
                        @test isequal(ds[name][:, 1], expected[name][:, 1])
                    end
                end
                @test ds["signal"][71, 1] == 3.0
                @test ds["signal_std"][71, 1] == 1.0
                @test ds["signal_error"][71, 1] == 10.0
                @test ds["signal_error_std"][71, 1] == 0.0
                @test ds["signal"].attrib["satellitegridding_std"] == "signal_std"
                @test ds["signal_std"].attrib["satellitegridding_mean"] == "signal"
            end
            tagged || continue
            Dataset(reordered, "a") do ds
                ds["signal_std"].attrib["satellitegridding_mean"] = "signal_error"
            end
            @test_throws ArgumentError coarsen_healpix(reordered, invalid)
            @test !ispath(invalid)
            Dataset(reordered, "a") do ds
                for name in ("n", "signal", "signal_error", "signal_std", "signal_error_std")
                    attrs = ds[name].attrib
                    for attr in ("satellitegridding_role", "satellitegridding_mean", "satellitegridding_std")
                        haskey(attrs, attr) && delete!(attrs, attr)
                    end
                end
            end
        end

        cfg.grid_vars["orphan_std"] = "error"
        grid_l2(cfg, spec, times; outfile=fine, compute_std=false, opts...)
        Dataset(fine, "a") do ds
            @test ds["signal_error"].attrib["comment"] == copied_comment
            delete!(ds.attrib, "accumulator_precision")
        end
        @test_logs (:warn, r"accumulator precision metadata is unavailable\. See the exact-archive") coarsen_healpix(fine, coarse)
        Dataset(fine, "a") do ds
            for name in ("n", "signal", "signal_error", "orphan_std")
                delete!(ds[name].attrib, "satellitegridding_role")
            end
            delete!(ds.attrib, "significant_digits_mean")
            delete!(ds.attrib, "geometry_precision")
        end
        @test_logs (:warn, r"rounding settings are unavailable") coarsen_healpix(fine, coarse)
        Dataset(coarse) do ds
            @test ds["signal"][71, 1] == 3.0
            @test ds["signal_error"][71, 1] == 10.0
            @test ds["signal_error"].attrib["satellitegridding_role"] == "mean"
            @test ds["orphan_std"][71, 1] == 10.0
            @test ds["orphan_std"].attrib["satellitegridding_role"] == "mean"
            @test !haskey(ds, "signal_std")
        end
        Dataset(fine, "a") do ds
            defVar(ds, "signal_std", Float64, ("time",))
            defVar(ds, "signal_error_std", Float64, ("cell", "time"))
        end
        @test_throws ArgumentError("Ambiguous untagged HEALPix layout: expected complete, disjoint <mean>_std pairs on (cell, time)") coarsen_healpix(fine, invalid; legacy_std=:by_name)
        @test !ispath(invalid)

        # Reproduce the legacy writer's mean-only layout without statistic tags.
        empty!(cfg.grid_vars)
        cfg.grid_vars["signal"] = "signal"
        cfg.grid_vars["signal_std"] = "error"
        grid_l2(cfg, spec, times; outfile=fine, compute_std=false, opts...)
        grid_l2(cfg, parent_spec, times; outfile=direct, compute_std=false, opts...)
        Dataset(fine, "a") do ds
            for name in ("n", "signal", "signal_std")
                delete!(ds[name].attrib, "satellitegridding_role")
            end
            for attr in ("accumulator_precision", "significant_digits_mean", "geometry_precision")
                delete!(ds.attrib, attr)
            end
        end
        ambiguity = r"Ambiguous untagged HEALPix archive.*legacy_std=:by_name.*legacy_std=:none"
        @test_throws ambiguity coarsen_healpix(fine, invalid)
        @test !ispath(invalid)
        args = ["coarsen", "--input", fine, "--output", invalid]
        @test_throws ambiguity ProgressTestCLI.main(args)
        @test !ispath(invalid)
        @test_throws ArgumentError("legacy_std must be :error, :by_name, or :none") coarsen_healpix(fine, invalid; legacy_std=:guess)
        @test !ispath(invalid)
        @test_logs (:warn, r"rounding settings are unavailable") ProgressTestCLI.main([args; "--legacy-std"; "none"])
        Dataset(invalid) do actual
            Dataset(direct) do expected
                for name in ("n", "signal", "signal_std")
                    @test isequal(actual[name][:, 1], expected[name][:, 1])
                end
            end
            @test actual["signal_std"][71, 1] == 10.0
            @test actual["signal_std"].attrib["satellitegridding_role"] == "mean"
            @test !haskey(actual["signal"].attrib, "satellitegridding_std")
        end
    end
end

@testset "Legacy zero-fill weights" begin
    mktempdir(joinpath(@__DIR__, "..", ".scratch")) do dir
        fine, coarse = joinpath.(dir, ("legacy.nc", "coarse.nc"))
        Dataset(fine, "c") do ds
            SatelliteGridding.write_coordinates!(ds, HEALPixGridSpec(nside=2))
            defDim(ds, "time", Inf)
            defVar(ds, "time", Float32, ("time",))[1] = 0f0
            signal = defVar(ds, "signal", Float32, ("cell", "time"); fillvalue=0f0)
            signal.var[1, 1] = 7f0
            deviation = defVar(ds, "signal_std", Float32, ("cell", "time"); fillvalue=-999f0)
            deviation.var[1, 1] = 0f0
            weights = defVar(ds, "n", Float32, ("cell", "time"); fillvalue=0f0)
            weights.var[:, 1] = [1f0; zeros(Float32, 47)]
        end
        @test_logs (:warn, r"rounding settings are unavailable") coarsen_healpix(fine, coarse; legacy_std=:by_name)
        Dataset(coarse) do ds
            @test ds["n"][:, 1] == [1.0; zeros(11)]
            @test ds["signal"][1, 1] == 7.0
            @test ds["signal_std"][1, 1] == 0.0
            @test all(ismissing, ds["signal"][2:end, 1])
        end
    end
end

@testset "Coarsened fill markers and footprint totals" begin
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        cfg = _healpix_ka_signal_config(dir, "fill")
        cfg.grid_vars["signal_std"] = "signal"  # A science field in a mean-only archive.
        cfg.options["output"] = Dict("precision" => "float64", "significant_digits_mean" => "none",
                                    "statistics" => ["sum_w2", "n_soundings"])
        lon1, lat1 = pix2ang_nest(8, 280)
        lon2, lat2 = pix2ang_nest(8, 281)
        whole(lon, lat) = (lon - 0.01, lon + 0.01, lat - 0.01, lat + 0.01)
        _write_healpix_ka_signal_input(joinpath(dir, "fill_2020-01-01.nc"),
            [whole(lon1, lat1), whole(lon2, lat2)], [-998.0, -1000.0])
        _write_healpix_ka_signal_input(joinpath(dir, "fill_2020-01-02.nc"),
            [(11.05, 11.45, -0.2, 0.2)], [3.0])
        times = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 2), Day(1))
        fine, coarse = joinpath.(dir, ("fine.nc", "coarse.nc"))
        grid_l2(cfg, HEALPixGridSpec(nside=8, T=Float64), times; outfile=fine,
                n_oversample=4, progress=false, verbose=false)
        coarsen_healpix(fine, coarse; levels=2, legacy_std=:by_name)
        Dataset(coarse) do ds
            parent = 18  # NESTED cell ID 17 at nside=2.
            @test ds["signal"][parent, 1] == ds["signal_std"][parent, 1] == -999.0
            @test ds["n"][parent, :] == [2.0, 1.0]
            @test ds["sum_w2"][parent, :] == [2.0, 5/16]
            @test ds["n_soundings"][parent, :] == [2, 4]
            @test !haskey(ds, "lon_bounds")
        end
        Dataset(fine, "a") do ds
            ds["n_soundings"].var[281, 1] = typemax(Int32)
        end
        @test_throws OverflowError("n_soundings exceeds the Int32 limit (2147483647)") coarsen_healpix(fine, coarse; levels=2)
    end
end
