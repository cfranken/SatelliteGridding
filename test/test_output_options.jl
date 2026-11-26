@testset "Output precision, rounding, and count fill" begin
    mktempdir() do dir
        # The first two footprints share a cell. The third occupies a different
        # cell. The fourth lies outside the rectangular grid and splits over
        # several HEALPix cells, including one with a single subpixel at n=3.
        # The binary values give exact Float64 means and standard deviations
        # and change when narrowed to Float32.
        mean_signal = 1 + 0x1.23456789abcdp-3
        deviation = 0x1.2345678p-5
        single_signal = 2 + 0x1.123456789abcp-2
        latitudes = [0.25, 0.25, 20.25, -40.0]
        longitudes = [0.25, 0.25, 20.25, 100.0]
        lat_halfwidth = [0.01, 0.01, 0.01, 18.0]
        lon_halfwidth = [0.01, 0.01, 0.01, 24.0]
        Dataset(joinpath(dir, "input_2020-01-01.nc"), "c") do ds
            defDim(ds, "obs", length(latitudes))
            defDim(ds, "vertex", 4)
            defVar(ds, "lat", Float64, ("obs",))[:] = latitudes
            defVar(ds, "lon", Float64, ("obs",))[:] = longitudes
            lat_bounds = defVar(ds, "lat_bounds", Float64, ("obs", "vertex"))
            lon_bounds = defVar(ds, "lon_bounds", Float64, ("obs", "vertex"))
            lat_bounds[:, :] = latitudes .+ lat_halfwidth .* [-1 -1 1 1]
            lon_bounds[:, :] = longitudes .+ lon_halfwidth .* [-1 1 1 -1]
            defVar(ds, "signal", Float64, ("obs",))[:] =
                [mean_signal - deviation, mean_signal + deviation, single_signal, 1.0]
        end
        config_file = joinpath(dir, "config.toml")
        write(config_file, """
        folder = "$dir"
        filePattern = "input_YYYY-MM-DD.nc"
        [basic]
        lat = "lat"
        lon = "lon"
        lat_bnd = "lat_bounds"
        lon_bnd = "lon_bounds"
        [grid]
        signal = "signal"
        [center]
        min_count = 1
        [output]
        precision = "float64"
        significant_digits_mean = "none"
        significant_digits_std = "none"
        significant_digits_center = "none"
        """)
        config = load_config(config_file)
        output = config.options["output"]
        rect = GridSpec(lat_min=0.0, lat_max=30.0, lon_min=0.0, lon_max=30.0, dlat=5.0, dlon=5.0)
        cube = CubedSphereGridSpec(Nc=8, T=Float64)
        healpix = HEALPixGridSpec(nside=4, T=Float64)
        times = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 2), Day(1))
        outfile = joinpath(dir, "output.nc")
        @test_throws ArgumentError grid_l2(config, rect, times; n_oversample=0, outfile)

        # A Float32 HEALPix grid with a Float64 KA accumulator must also retain
        # the full accumulator precision when Float64 output is requested.
        for (spec, backend) in ((rect, nothing), (rect, CPU()), (cube, nothing),
                                (healpix, nothing), (healpix, CPU()),
                                (HEALPixGridSpec(nside=4, T=Float32), CPU()))
            compute_std = !(spec === rect && backend !== nothing)
            output["fill_value"] = -1.234567890123e100
            grid_l2(config, spec, times; backend, compute_std, outfile, n_oversample=3)
            Dataset(outfile) do ds
                for key in (compute_std ? ("signal", "signal_std", "n") : ("signal", "n"))
                    @test eltype(ds[key].var) === Float64
                    @test NCDatasets.deflate(ds[key].var) == (false, true, 4)
                    expected_fill = key == "n" ? -1.0 : output["fill_value"]
                    @test ds[key].attrib["_FillValue"] === expected_fill
                end
                weights = ds["n"][:]
                if spec isa HEALPixGridSpec && backend !== nothing
                    @test any(==(1 / 9), weights)
                end
                @test count(==(1.0), weights) == count(==(2.0), weights) == 1
                @test ds["signal"][:][weights .== 1.0] == [single_signal]
                @test ds["signal"][:][weights .== 2.0] == [mean_signal]
                @test all(ismissing, ds["signal"][:][iszero.(weights)])
                if compute_std
                    @test ds["signal_std"][:][weights .== 1.0] == [0.0]
                    @test ds["signal_std"][:][weights .== 2.0] == [deviation]
                end
            end
        end

        output["precision"] = "float32"
        output["significant_digits_mean"] = 4
        output["significant_digits_std"] = 3
        for fill in (0, 1), spec in (rect, cube, healpix)
            output["fill_value"] = fill
            grid_l2(config, spec, times; compute_std=true, outfile, n_oversample=4)
            Dataset(outfile) do ds
                @test ds["n"].attrib["_FillValue"] === -1f0
                @test ds["signal"].attrib["_FillValue"] === Float32(fill)
                weights = ds["n"][:]
                @test !any(ismissing, weights)
                @test any(iszero, weights)
                @test count(==(1f0), weights) == 1
                @test ds["signal"][:][weights .== 1f0] == [round(Float32(single_signal); sigdigits=4)]
                @test ds["signal_std"][:][weights .== 2f0] == [round(Float32(deviation); sigdigits=3)]
            end
        end

        output["precision"] = "float64"
        output["significant_digits_mean"] = "none"
        output["significant_digits_std"] = "none"
        rolling = TimeSpec(times.start_date, times.stop_date, Day(1);
                           sample_step=Day(1), window_halfwidth_days=1)
        grid_l2(config, HEALPixGridSpec(nside=4, T=Float32), rolling;
                backend=CPU(), compute_std=true, outfile, n_oversample=3)
        Dataset(outfile) do ds
            @test all(slice -> any(==(1 / 9), slice), eachcol(ds["n"][:, :]))
            weights = ds["n"][:]
            @test ds["signal"][:][weights .== 2.0] == fill(mean_signal, 2)
            @test ds["signal_std"][:][weights .== 2.0] == fill(deviation, 2)
        end

        grid_center(config, rect, times; outfile, geo_provider=:variables)
        Dataset(outfile) do ds
            @test eltype(ds["signal"].var) === Float64
            @test ds["n"].attrib["_FillValue"] === -1.0
            weights = ds["n"][:]
            @test !any(ismissing, weights)
            @test count(==(1.0), weights) == 1
            @test ds["signal"][:][weights .== 2.0] == [mean_signal]
        end

        nc_path = joinpath(dir, "rejected.nc")
        zarr_path = joinpath(dir, "rejected.zarr")
        @test_throws "Zarr output requires output.precision" grid_l2(
            config, rect, times; outfile=nc_path, zarr_out=zarr_path)
        @test !ispath(nc_path) && !ispath(zarr_path)
    end
end

@testset "Output rounding at precision boundaries" begin
    cases = (
        (floatmax(Float32), "float64", 4, 3.403e38),
        (1.234375f0, "float64", 3, 1.23),
        (1.75e308, "float64", 1, missing),
        (floatmax(Float32), "float32", 4, missing),
        (2 * Float64(floatmax(Float32)), "float32", "none", missing),
        (1.23499999999, "float32", 3, 1.24f0),
        (1.0737439f0, "float32", 10, 1.0737439f0),
        (π / 4, "float64", 20, π / 4),
        (0.001, "float64", typemax(Int), 0.001),
    )
    mktempdir() do dir
        grid_vars = OrderedDict("signal" => "signal")
        config = DataSourceConfig(Dict{String,String}(), grid_vars, FilterRule[], "", dir)
        spec = GridSpec(lat_min=0f0, lat_max=1f0, lon_min=0f0, lon_max=1f0)
        outfile = joinpath(dir, "rounding.nc")
        date = DateTime(2020, 1, 1)
        for (value, precision, digits, expected) in cases
            options = NetCDFOptions(; precision, significant_digits_mean=digits,
                                     significant_digits_std=digits, significant_digits_center=digits)
            weights = ones(typeof(value), 1, 1)
            data = fill(value, 1, 1, 1)
            ds, vars = create_output_dataset(outfile, spec, 2, grid_vars, true; output_options=options)
            try
                SatelliteGridding._write_time_slice!(ds, vars, config, spec, options,
                    data, data, weights, true, 1, date)
                SatelliteGridding._write_center_time_slice!(ds, vars, config, options, data,
                    weights, false, 2, date + Day(1), 1)
            finally
                close(ds)
            end
            Dataset(outfile) do ds
                @test ds["signal"][1, 1, 1] === expected
                @test ds["signal_std"][1, 1, 1] === expected
                @test ds["signal"][1, 1, 2] === expected
                @test ds["n"][1, 1, :] == [1, 1]
            end
        end
    end

    # Narrowing before the legacy Float32 mean rounding changes this value by one ulp.
    data = [1.0747630006056146]
    digits = NetCDFOptions().significant_digits_mean
    @test SatelliteGridding._masked_output_slice(data, [1.0], -999f0, digits, 8) == [1.074763f0]
end
