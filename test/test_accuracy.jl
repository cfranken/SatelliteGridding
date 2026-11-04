@testset "Accuracy — analytic synthetic inputs through NetCDF" begin
    @testset "Four observations retain their own corners in either bounds layout" begin
        mktempdir() do dir
            rectangles = [(x - 0.5, x + 0.5, y - 0.5, y + 0.5)
                          for (x, y) in zip((-100, -50, 50, 100), (-30, -10, 10, 30))]
            signals = Float32[10, 20, 30, 40]
            gs = GridSpec(lat_min=-40f0, lat_max=40f0, lon_min=-120f0, lon_max=120f0,
                          dlat=20f0, dlon=60f0)
            ts = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 1), Dates.Day(1))
            input = joinpath(dir, "four_2020-01-01.nc")
            expected_weights = zeros(Float32, 4, 4)
            for i in 1:4
                expected_weights[i, i] = 1
            end
            for vertex_first in (false, true), grouped in (false, true)
                _accuracy_input(input, rectangles, signals; vertex_first=vertex_first,
                                bounds_group=grouped ? "Geolocation" : nothing)
                prefix = grouped ? "Geolocation/" : ""
                cfg = DataSourceConfig(Dict("lat_bnd" => prefix * "lat_bounds",
                                            "lon_bnd" => prefix * "lon_bounds"),
                    OrderedDict("signal" => "signal"), FilterRule[], "four_YYYY-MM-DD.nc", dir)
                for backend in (nothing, CPU())
                    out = joinpath(dir, "four_output.nc")
                    grid_l2(cfg, gs, ts; backend=backend, outfile=out, n_oversample=4)
                    Dataset(out) do ds
                        @test ds["n"][:, :, 1] ≈ expected_weights
                        for i in 1:4
                            @test ds["signal"][i, i, 1] ≈ signals[i]
                        end
                        @test count(ismissing, ds["signal"][:, :, 1]) == 12
                    end
                end
                out = joinpath(dir, "four_cube.nc")
                grid_l2(cfg, CubedSphereGridSpec(Nc=8), ts; outfile=out, n_oversample=4)
                Dataset(out) do ds
                    @test sum(ds["n"][:, :, :, 1]) ≈ 4
                    @test Set(skipmissing(vec(ds["signal"][:, :, :, 1]))) == Set(signals)
                end
            end
            # Selecting rows must retain the previous rejection of science
            # variables whose full length does not match the sounding axis.
            for n_bad in (2, 5)
                path = "bad_signal_$n_bad"
                Dataset(input, "a") do ds
                    defDim(ds, "bad_dim_$n_bad", n_bad)
                    defVar(ds, path, Float32, ("bad_dim_$n_bad",))[:] = ones(Float32, n_bad)
                end
                cfg = DataSourceConfig(Dict("lat_bnd" => "Geolocation/lat_bounds",
                                            "lon_bnd" => "Geolocation/lon_bounds"),
                    OrderedDict("signal" => path), FilterRule[], "four_YYYY-MM-DD.nc", dir)
                @test_throws ErrorException grid_l2(cfg, gs, ts;
                    outfile=joinpath(dir, "invalid.nc"), n_oversample=4)
            end
        end
    end

    mktempdir() do dir
        broad = (0.5, 1.5, 0.25, 0.75) # half its area in each of two cells
        left = (0.1, 0.9, 0.25, 0.75)  # entirely in the left cell
        _accuracy_input(joinpath(dir, "fake_2020-01-01.nc"),
                        [broad, left, left, left, left],
                        Float32[10, 1000, -9999, NaN, Inf]; qa=Float32[1, 0, 1, 1, 1], version="B1")
        _accuracy_input(joinpath(dir, "fake_2020-01-02.nc"), [left], Float32[30]; version="B2")
        config_file = joinpath(dir, "config.toml")
        write(config_file, """
        filePattern = "fake_YYYY-MM-DD.nc"
        folder = "$(replace(dir, '\\' => '/'))"
        [basic]
        lat = "lat"
        lon = "lon"
        lat_bnd = "lat_bounds"
        lon_bnd = "lon_bounds"
        [grid]
        signal = "signal"
        packed = "packed"
        [filter]
        qa = "> 0.5"
        [metadata]
        attributes = ["@institution", "Metadata/BuildId", "@absent"]
        variable_attributes = ["units", "long_name"]
        [output]
        shuffle = true
        fill_value = nan
        chunksizes = [2, 1, 1]
        checksum = "fletcher32"
        """)
        config = load_config(config_file)
        gs = GridSpec(lat_min=0f0, lat_max=2f0, lon_min=0f0, lon_max=3f0, dlat=1f0, dlon=1f0)
        ts = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 3), Dates.Day(2))

        # Sequential and KA CPU must agree with the analytic result, not merely
        # with each other. Two snapshots: both days combined, then no input.
        for backend in (nothing, CPU())
            out = joinpath(dir, backend === nothing ? "sequential.nc" : "ka.nc")
            grid_l2(config, gs, ts; outfile=out, backend=backend,
                    n_oversample=40, compute_std=backend === nothing)
            Dataset(out) do ds
                w = ds["n"][:, :, 1]
                @test w ≈ Float32[1.5 0; 0.5 0; 0 0] atol=3f-5
                @test sum(w) ≈ 2f0 atol=3f-5
                @test ds["signal"][1, 1, 1] ≈ 70/3 atol=2e-4
                @test ds["signal"][2, 1, 1] ≈ 10 atol=2e-4
                @test ds["packed"][1, 1, 1] ≈ 161/3 atol=4e-4
                @test ds["packed"][2, 1, 1] ≈ 27 atol=4e-4
                if backend === nothing
                    @test ds["signal_std"][1, 1, 1] ≈ sqrt(800/9) atol=2e-4
                    @test ds["packed_std"][1, 1, 1] ≈ 2sqrt(800/9) atol=3e-4
                    @test ds["signal_std"][2, 1, 1] ≈ 0 atol=1e-6
                    @test all(ismissing, ds["signal_std"][:, :, 2])
                end
                @test ismissing(ds["signal"][3, 1, 1])
                @test all(ismissing, ds["signal"][:, :, 2])
                @test all(iszero, ds["n"][:, :, 2])
                @test isnan(ds["signal"].attrib["_FillValue"])
                @test !haskey(ds["packed"].attrib, "scale_factor")
                @test ds["signal"].attrib["units"] == "ppm"
                @test ds.attrib["satellite_product_version"] == "B1; B2"
                versions = SatelliteGridding.JSON.parse(ds.attrib["source_versions"])
                @test versions["Metadata/BuildId"] == ["B1.build", "B2.build"]
                attrs = SatelliteGridding.JSON.parse(ds.attrib["source_metadata"])
                @test attrs["@institution"] == ["Synthetic lab"]
                @test !haskey(attrs, "@absent")
                @test ds["time"][:] == [DateTime(2020, 1, 1), DateTime(2020, 1, 3)]
                @test NCDatasets.deflate(ds["signal"].var) == (true, true, 4)
                @test NCDatasets.chunking(ds["signal"].var)[2] == (2, 1, 1)
            end
        end

        # Rolling windows reuse daily accumulations. Day 3 includes only day 2;
        # the right cell must become missing rather than keep day 1's value.
        rolling = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 3), Dates.Day(1);
                           sample_step=Dates.Day(1), window_halfwidth_days=1)
        out = joinpath(dir, "rolling.nc")
        grid_l2(config, gs, rolling; outfile=out, n_oversample=40, compute_std=true)
        Dataset(out) do ds
            @test ds["signal"][1, 1, :] ≈ [70/3, 70/3, 30] atol=2e-4
            @test ds["n"][1, 1, :] ≈ [1.5, 1.5, 1] atol=3e-5
            @test ismissing(ds["signal"][2, 1, 3])
            @test ds["signal_std"][1, 1, 3] ≈ 0 atol=1e-5
            @test ds.attrib["satellite_product_version"] == "B1; B2"
        end

        # Center mode: arithmetic mean after unpacking, min_count, and NaN fill.
        _accuracy_input(joinpath(dir, "center_2020-01-01.nc"),
                        [left, left, left, left], Float32[10, 30, 1000, -9999];
                        qa=Float32[1, 1, 0, 1], version="B3")
        config.options["min_count"] = 2
        center_config = DataSourceConfig(config.basic, config.grid_vars, config.filters,
                                         "center_YYYY-MM-DD.nc", dir, config.options)
        out = joinpath(dir, "center.nc")
        grid_center(center_config, gs, ts; outfile=out, geo_provider=:variables)
        Dataset(out) do ds
            @test ds["signal"][1, 1, 1] ≈ 20
            @test ds["packed"][1, 1, 1] ≈ 47
            @test ds["n"][1, 1, 1] == 2
            @test ismissing(ds["signal"][2, 1, 1])
            @test all(ismissing, ds["signal"][:, :, 2])
            @test ds.attrib["satellite_product_version"] == "B3"
        end
    end
end
