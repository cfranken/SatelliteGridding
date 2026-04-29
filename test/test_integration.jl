@testset "Integration — end-to-end L2 gridding" begin

    @testset "Synthetic L2 file → grid_l2 pipeline" begin
        test_dir = mktempdir()

        # --- Create synthetic L2 NetCDF file ---
        # Small 10°×10° domain, 2° resolution → 5×5 grid
        # 20 footprints with known positions and signal values
        n_soundings = 20
        input_file = joinpath(test_dir, "TROPO_SIF_2020-01-01_test.nc")

        ds = Dataset(input_file, "c")
        defDim(ds, "x", n_soundings)
        defDim(ds, "nv", 4)

        # Center coordinates: scattered across the domain
        lats_center = Float32[1, 2, 3, 4, 5, 6, 7, 8, 9, 1,
                              2, 3, 4, 5, 6, 7, 8, 9, 5, 5]
        lons_center = Float32[1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
                              8, 7, 6, 5, 4, 3, 2, 1, 5, 5]

        # Create small footprints (0.5° × 0.5°) around each center
        # All will fit within a single 2° grid cell (fast path)
        lat_bnds = zeros(Float32, n_soundings, 4)
        lon_bnds = zeros(Float32, n_soundings, 4)
        for i in 1:n_soundings
            dx = 0.25f0
            # Corners: SW, SE, NE, NW
            lat_bnds[i, :] = [lats_center[i] - dx, lats_center[i] - dx,
                              lats_center[i] + dx, lats_center[i] + dx]
            lon_bnds[i, :] = [lons_center[i] - dx, lons_center[i] + dx,
                              lons_center[i] + dx, lons_center[i] - dx]
        end

        # Signal = lat + lon (known analytical function)
        signal = lats_center .+ lons_center

        v_lat = defVar(ds, "lat", Float32, ("x",))
        v_lon = defVar(ds, "lon", Float32, ("x",))
        v_lat_bnd = defVar(ds, "lat_bnds", Float32, ("x", "nv"))
        v_lon_bnd = defVar(ds, "lon_bnds", Float32, ("x", "nv"))
        v_sif = defVar(ds, "sif", Float32, ("x",),
                       attrib=["units" => "mW/m²/sr/nm", "long_name" => "SIF"])
        v_sza = defVar(ds, "sza", Float32, ("x",))

        v_lat[:] = lats_center
        v_lon[:] = lons_center
        v_lat_bnd[:, :] = lat_bnds
        v_lon_bnd[:, :] = lon_bnds
        v_sif[:] = signal
        v_sza[:] = fill(Float32(30.0), n_soundings)  # all pass sza < 80 filter

        close(ds)

        # --- Create TOML config ---
        config_file = joinpath(test_dir, "test_config.toml")
        config_toml = """
        filePattern = "TROPO_SIF_YYYY-MM-DD_*.nc"
        folder = "$(replace(test_dir, "\\" => "/"))/"

        [basic]
        lat = "lat"
        lon = "lon"
        lat_bnd = "lat_bnds"
        lon_bnd = "lon_bnds"

        [grid]
        sif = "sif"

        [filter]
        sza = "< 80"
        """
        write(config_file, config_toml)

        # --- Run through dispatch API ---
        output_file = joinpath(test_dir, "output.nc")
        config = load_config(config_file)
        grid_spec = GridSpec(lat_min=0.0f0, lat_max=10.0f0,
                             lon_min=0.0f0, lon_max=10.0f0,
                             dlat=2.0f0, dlon=2.0f0)
        time_spec = TimeSpec(DateTime("2020-01-01"), DateTime("2020-01-01"),
                             Dates.Day(1); oversample_temporal=1.0f0)

        grid(config, grid_spec, time_spec, SubpixelGridding(n_oversample=10);
             compute_std=true, outfile=output_file)

        # --- Verify output ---
        ds_out = Dataset(output_file)

        @test size(ds_out["lat"]) == (5,)
        @test size(ds_out["lon"]) == (5,)
        @test size(ds_out["time"]) == (1,)

        weights = ds_out["n"][1, :, :]
        sif_grid = ds_out["sif"][1, :, :]

        # At least some cells should have data
        @test maximum(weights) > 0

        # Cells with data should have reasonable SIF values (between 2 and 18)
        for i in eachindex(sif_grid)
            if weights[i] > 0 && sif_grid[i] > -999
                @test 1.0f0 < sif_grid[i] < 19.0f0
            end
        end

        # The center cell (lon=4-6, lat=4-6, center at 5,5) has footprints 19 and 20
        # both with signal = 5+5 = 10. Check that cell.
        # Grid indices: lon index 3 (covers 4-6°), lat index 3 (covers 4-6°)
        center_lon_idx = 3
        center_lat_idx = 3
        if weights[center_lon_idx, center_lat_idx] > 0
            @test sif_grid[center_lon_idx, center_lat_idx] ≈ 10.0f0 rtol=0.1
        end

        # Check that std was computed (should exist and be non-negative where data exists)
        sif_std = ds_out["sif_std"][1, :, :]
        for i in eachindex(sif_std)
            if weights[i] > 0 && sif_std[i] > -999
                @test sif_std[i] >= 0
            end
        end

        close(ds_out)

        # Cleanup
        rm(test_dir; recursive=true)
    end

    @testset "Synthetic center-coordinate file → grid_center pipeline" begin
        test_dir = mktempdir()
        input_file = joinpath(test_dir, "CENTER_2020-01-01.nc")

        ds = Dataset(input_file, "c")
        defDim(ds, "y", 2)
        defDim(ds, "x", 2)
        v_lat = defVar(ds, "lat", Float32, ("y", "x"))
        v_lon = defVar(ds, "lon", Float32, ("y", "x"))
        v_band = defVar(ds, "band", Float32, ("y", "x"))

        v_lat[:, :] = Float32[0.5 0.5; 2.5 2.5]
        v_lon[:, :] = Float32[0.5 2.5; 0.5 2.5]
        v_band[:, :] = Float32[10 20; 30 40]
        close(ds)

        config_file = joinpath(test_dir, "center_config.toml")
        config_toml = """
        filePattern = "CENTER_YYYY-MM-DD.nc"
        folder = "$(replace(test_dir, "\\" => "/"))/"

        [basic]
        lat = "lat"
        lon = "lon"

        [center]
        min_count = 1

        [grid]
        band = "band"
        """
        write(config_file, config_toml)

        output_file = joinpath(test_dir, "center_output.nc")
        config = load_config(config_file)
        grid_spec = GridSpec(lat_min=0.0f0, lat_max=4.0f0,
                             lon_min=0.0f0, lon_max=4.0f0,
                             dlat=2.0f0, dlon=2.0f0)
        time_spec = TimeSpec(DateTime("2020-01-01"), DateTime("2020-01-01"),
                             Dates.Day(1))

        grid(config, grid_spec, time_spec, CenterPointGridding(); outfile=output_file)

        ds_out = Dataset(output_file)
        weights = ds_out["n"][1, :, :]
        band = ds_out["band"][1, :, :]

        @test weights ≈ Float32[1 1; 1 1]
        @test band[1, 1] ≈ 10.0f0
        @test band[2, 1] ≈ 20.0f0
        @test band[1, 2] ≈ 30.0f0
        @test band[2, 2] ≈ 40.0f0

        close(ds_out)
        rm(test_dir; recursive=true)
    end

    @testset "Generated MODIS geolocation cache → grid_center pipeline" begin
        test_dir = mktempdir()
        input_file = joinpath(test_dir, "MCD43A4.A2020001.h18v09.test.nc")

        lat_ul, lon_ul = modis_sinusoidal_latlon(8, 4, 1, 1)
        @test lat_ul ≈ 49.997917f0 atol=1f-5
        @test lon_ul ≈ -155.5624f0 atol=1f-4

        ds = Dataset(input_file, "c")
        defDim(ds, "y", 2)
        defDim(ds, "x", 2)
        v_band = defVar(ds, "Nadir_Reflectance_Band1", Float32, ("y", "x"))
        v_band[:, :] = Float32[10 20; 30 40]
        close(ds)

        config_file = joinpath(test_dir, "modis_center_config.toml")
        config_toml = """
        filePattern = "MCD43A4.AYYYYDOY.*.nc"
        folder = "$(replace(test_dir, "\\" => "/"))/"

        [center]
        min_count = 1
        modis_pixels = 2

        [grid]
        band = "Nadir_Reflectance_Band1"
        """
        write(config_file, config_toml)

        output_file = joinpath(test_dir, "modis_center_output.nc")
        geo_cache = joinpath(test_dir, "geo_cache")
        config = load_config(config_file)
        grid_spec = GridSpec(lat_min=-10.0f0, lat_max=0.0f0,
                             lon_min=0.0f0, lon_max=10.0f0,
                             dlat=5.0f0, dlon=5.0f0)
        time_spec = TimeSpec(DateTime("2020-01-01"), DateTime("2020-01-01"),
                             Dates.Day(1))

        grid(config, grid_spec, time_spec, CenterPointGridding();
             geo_cache=geo_cache,
             geo_provider=:modis,
             outfile=output_file)

        @test isfile(modis_tile_cache_path(geo_cache, 18, 9; pixels=2))

        ds_out = Dataset(output_file)
        weights = ds_out["n"][1, :, :]
        band = ds_out["band"][1, :, :]

        @test weights ≈ Float32[1 1; 1 1]
        @test band[1, 1] ≈ 30.0f0
        @test band[2, 1] ≈ 40.0f0
        @test band[1, 2] ≈ 10.0f0
        @test band[2, 2] ≈ 20.0f0

        close(ds_out)
        rm(test_dir; recursive=true)
    end

    @testset "File pattern resolution" begin
        @test SatelliteGridding._resolve_date_placeholders(
            "data_YYYY-MM-DD.nc", DateTime("2020-03-15")) == "data_2020-03-15.nc"

        @test SatelliteGridding._resolve_date_placeholders(
            "YYYY/DOY/file.nc", DateTime("2020-02-01")) == "2020/032/file.nc"
    end
end
