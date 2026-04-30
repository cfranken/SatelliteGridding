@testset "Configuration" begin

    @testset "load_config — TOML with string filter expressions" begin
        toml_content = """
        filePattern = "data_YYYY-MM-DD_*.nc"
        folder = "/tmp/test_data/"

        [basic]
        lat = "latitude"
        lon = "longitude"
        lat_bnd = "lat_corners"
        lon_bnd = "lon_corners"

        [grid]
        sif = "SIF_757nm"
        reflectance = "continuum_radiance"

        [filter]
        qa_value = "> 0.5"
        sza = "< 80"
        mode = "== 1"
        ch4 = "1600 < x < 2200"
        """
        tmpfile = tempname() * ".toml"
        write(tmpfile, toml_content)

        config = load_config(tmpfile)

        @test config.basic["lat"] == "latitude"
        @test config.basic["lon_bnd"] == "lon_corners"
        @test length(config.grid_vars) == 2
        @test config.grid_vars["sif"] == "SIF_757nm"
        @test config.file_pattern == "data_YYYY-MM-DD_*.nc"
        @test config.folder == "/tmp/test_data/"

        # Check filter rules
        @test length(config.filters) == 4
        rules = Dict(r.variable => r for r in config.filters)
        @test rules["qa_value"].op == :gt
        @test rules["qa_value"].lo == 0.5
        @test rules["sza"].op == :lt
        @test rules["sza"].lo == 80.0
        @test rules["mode"].op == :eq
        @test rules["mode"].lo == 1.0
        @test rules["ch4"].op == :between
        @test rules["ch4"].lo == 1600.0
        @test rules["ch4"].hi == 2200.0

        rm(tmpfile)
    end

    @testset "load_config — JSON backward compatibility (legacy filters)" begin
        json_content = """
        {
            "basic": {"lat": "lat", "lon": "lon", "lat_bnd": "lat_b", "lon_bnd": "lon_b"},
            "grid": {"var1": "path/to/var1"},
            "filter_gt": {"qa": 0.5},
            "filter_lt": {"sza": 80},
            "filter_eq": {"mode": 1},
            "filePattern": "*.nc",
            "folder": "/data/"
        }
        """
        tmpfile = tempname() * ".json"
        write(tmpfile, json_content)

        config = load_config(tmpfile)

        @test length(config.filters) == 3
        rules = Dict(r.variable => r for r in config.filters)
        @test rules["qa"].op == :gt
        @test rules["qa"].lo == 0.5
        @test rules["sza"].op == :lt
        @test rules["sza"].lo == 80.0
        @test rules["mode"].op == :eq
        @test rules["mode"].lo == 1.0

        rm(tmpfile)
    end

    @testset "load_config — minimal TOML (no filters)" begin
        toml_content = """
        filePattern = "*.nc"
        folder = "/data/"

        [basic]
        lat = "lat"
        lon = "lon"
        lat_bnd = "lat_b"
        lon_bnd = "lon_b"

        [grid]
        var1 = "path/to/var1"
        """
        tmpfile = tempname() * ".toml"
        write(tmpfile, toml_content)

        config = load_config(tmpfile)

        @test isempty(config.filters)
        @test length(config.grid_vars) == 1

        rm(tmpfile)
    end

    @testset "load_config — center config without basic section" begin
        toml_content = """
        filePattern = "*.hdf"
        folder = "/data/"

        [center]
        scale_factor = 0.0001

        [grid]
        var1 = "Nadir_Reflectance_Band1"
        """
        tmpfile = tempname() * ".toml"
        write(tmpfile, toml_content)

        config = load_config(tmpfile)

        @test isempty(config.basic)
        @test config.options["scale_factor"] == 0.0001
        @test config.grid_vars["var1"] == "Nadir_Reflectance_Band1"

        rm(tmpfile)
    end

    @testset "load_config — circular radius options" begin
        toml_content = """
        filePattern = "*.nc"
        folder = "/data/"

        [basic]
        lat = "Latitude"
        lon = "Longitude"

        [circle]
        radius = 5.25
        radius_unit = "km"

        [grid]
        sif = "Daily_Averaged_SIF"
        """
        tmpfile = tempname() * ".toml"
        write(tmpfile, toml_content)

        config = load_config(tmpfile)

        @test config.basic["lat"] == "Latitude"
        @test config.options["radius"] == 5.25
        @test config.options["radius_unit"] == "km"

        rm(tmpfile)
    end

    @testset "load_config — missing required section" begin
        toml_content = """
        [basic]
        lat = "lat"
        """
        tmpfile = tempname() * ".toml"
        write(tmpfile, toml_content)

        @test_throws ErrorException load_config(tmpfile)

        rm(tmpfile)
    end

    @testset "filter expression parsing" begin
        parse = SatelliteGridding._parse_filter_expr

        r = parse("sza", "< 80")
        @test r.op == :lt && r.lo == 80.0

        r = parse("qa", "> 0.5")
        @test r.op == :gt && r.lo == 0.5

        r = parse("mode", "== 1")
        @test r.op == :eq && r.lo == 1.0

        r = parse("ch4", "1600 < x < 2200")
        @test r.op == :between && r.lo == 1600.0 && r.hi == 2200.0

        r = parse("ch4", "1.6e3 < x < 2.2e3")
        @test r.op == :between && r.lo ≈ 1600.0 && r.hi ≈ 2200.0

        @test_throws ErrorException parse("var", "not a filter")
    end

    @testset "GridSpec construction" begin
        gs = GridSpec(lat_min=-10.0f0, lat_max=10.0f0,
                      lon_min=-20.0f0, lon_max=20.0f0,
                      dlat=5.0f0, dlon=10.0f0)

        @test gs.lat ≈ Float32[-7.5, -2.5, 2.5, 7.5]
        @test gs.lon ≈ Float32[-15.0, -5.0, 5.0, 15.0]
        @test length(gs.lat) == 4
        @test length(gs.lon) == 4
    end

    @testset "TimeSpec construction" begin
        ts = TimeSpec(DateTime("2020-01-01"), DateTime("2020-12-31"),
                      Dates.Day(8); oversample_temporal=2.0f0)

        @test ts.start_date == DateTime("2020-01-01")
        @test ts.oversample_temporal == 2.0f0
        @test ts.time_step == Dates.Day(8)
    end

    @testset "CLI argument helpers" begin
        args = Dict(
            "latMin" => -45.0f0, "latMax" => 45.0f0,
            "lonMin" => -90.0f0, "lonMax" => 90.0f0,
            "dLat" => 2.0f0, "dLon" => 2.0f0,
            "startDate" => "2020-06-01", "stopDate" => "2020-06-30",
            "dDays" => 8, "monthly" => false,
            "oversample_temporal" => 1.0f0,
        )

        gs = args_to_grid_spec(args)
        @test length(gs.lat) == 45
        @test length(gs.lon) == 90

        ts = args_to_time_spec(args)
        @test ts.start_date == DateTime("2020-06-01")
        @test ts.time_step == Dates.Day(8)
    end
end
