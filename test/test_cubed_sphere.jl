using SHA

function _cubed_sphere_accumulator_bits(spec::CubedSphereGridSpec{T}) where {T}
    # Repeated, overlapping footprints across a panel edge and the dateline,
    # plus a near-pole footprint and one wholly inside a cell.
    lat = [-2.2 -2.1 1.9 2.6; -1.4 -1.3 2.0 2.2; 3.0 3.1 7.6 7.0;
           4.0 4.0 8.0 8.0; 86.0 87.0 89.0 89.2; -0.02 -0.02 0.02 0.02]
    lon = [32.1 49.3 47.7 34.0; 33.2 50.0 46.8 37.0; 176.1 -176.2 -177.0 177.4;
           177.2 -178.0 -177.3 177.4; 127.0 -137.0 -101.0 113.0;
           -10.02 -9.98 -9.98 -10.02]
    values = T[1.1 9.2; -4.3 6.4; 13.5 -2.6; 31.7 5.8; 0.9 4.1; 8.2 1.3]
    n = 7
    shape = accumulator_shape(spec)
    data, m2, weights = zeros(T, shape..., 2), zeros(T, shape..., 2), zeros(T, shape)
    observations = (lat_bounds=lat, lon_bounds=lon, values=values)
    SatelliteGridding._accumulate_l2_spherical!(SatelliteGridding._spherical_subdivision(spec),
        data, m2, weights, true, spec, observations, n)
    # Serialize the raw bits in a fixed byte order before hashing.
    bits_type = T == Float32 ? UInt32 : UInt64
    bits = reinterpret(bits_type, vcat(vec(data), vec(m2), vec(weights)))
    bytes2hex(sha256(reinterpret(UInt8, htol.(bits))))
end
@testset "Cubed-sphere grid" begin


    @testset "Bit-identical sequential accumulation" begin
        # SHA-256 of the raw mean, M2, and weight bits from the pre-refactor
        # accumulators at 2c5e505. The 7×7 samples use non-binary weights.
        expected = (
            ("a88db4aa8b05ae83f8e0de9ced262d1eeeb4fb0a63a9854ba0be948dca89a497",
             "848da16b7cfee56a4425e5db059e481c4e018cdc6b6f606d35388f01f6999b2d"),
            ("7332b6b4bdb47a3f6ed357d7e4e1917110e1f38d321988d38be6696a87f5897a",
             "a359297feb40f50793aba8f034d14f76c3fc0c652204e60c564536f9339ec896"),
        )
        definitions = (GMAOCubedSphereDefinition(), EquiangularCubedSphereDefinition())
        for (d, def) in enumerate(definitions), (t, T) in enumerate((Float32, Float64))
            spec = CubedSphereGridSpec(Nc=9, definition=def, T=T)
            @test SatelliteGridding._spherical_subpixel_type(spec) === T
            @test _cubed_sphere_accumulator_bits(spec) == expected[d][t]
        end
    end

    @testset "Projection round-trip + whole-sphere coverage" begin
        for def in (GMAOCubedSphereDefinition(), EquiangularCubedSphereDefinition())
            spec = CubedSphereGridSpec(Nc=24, definition=def, T=Float64)
            Nc = spec.Nc

            @test grid_shape(spec) == (Nc, 6Nc)
            @test accumulator_shape(spec) == (Nc, 6Nc)
            @test size(spec.centers_lon) == (Nc, Nc, 6)
            @test size(spec.corners_lon) == (Nc + 1, Nc + 1, 6)

            # Every cell center inverts back to its own (i, j, panel).
            mismatches = 0
            for p in 1:6, j in 1:Nc, i in 1:Nc
                lon = spec.centers_lon[i, j, p]
                lat = spec.centers_lat[i, j, p]
                panel, s, t = lonlat_to_panel_xy(spec.definition, Nc, lon, lat, Float64)
                (panel == p && floor(Int, s) == i && floor(Int, t) == j) || (mismatches += 1)
            end
            @test mismatches == 0

            # Dense lon/lat sample: every point lands in exactly one valid folded cell.
            out_of_range = 0
            for lat in -89.0:1.7:89.0, lon in 0.0:1.9:359.0
                col, row = SatelliteGridding.to_fractional_index(spec, lat, lon)
                c = floor(Int, col)
                r = floor(Int, row)
                (1 <= c <= Nc && 1 <= r <= 6Nc) || (out_of_range += 1)
            end
            @test out_of_range == 0
        end
    end

    @testset "GMAO definition matches the GEOS/GCHP convention" begin
        spec = CubedSphereGridSpec(Nc=48)  # defaults to GMAOCubedSphereDefinition
        @test spec.definition isa CubedSphereDefinition
        @test SatelliteGridding.longitude_offset_deg(spec.definition) == -10.0
        @test SatelliteGridding.coordinate_law(spec.definition) isa GMAOEqualDistanceGnomonic
        @test SatelliteGridding.center_law(spec.definition) isa FourCornerNormalizedCenter
        @test eltype(spec.centers_lon) == Float32
        # Longitudes in [0, 360); latitudes in [-90, 90].
        @test all(0 .<= spec.centers_lon .< 360.0001)
        @test all(-90.0001 .<= spec.centers_lat .<= 90.0001)
    end

    @testset "Synthetic L2 → grid_l2 cubed-sphere pipeline" begin
        test_dir = mktempdir()
        spec = CubedSphereGridSpec(Nc=6, T=Float32)
        Nc = spec.Nc

        # One small footprint per panel, centered on cell (3, 3) of that panel,
        # with a distinct signal value. Small (±0.05°) so each is sub-cell
        # (fast path) at C6 (~15° cells).
        cells = [(3, 3, p) for p in 1:6]
        n = length(cells)
        lats_center = Float32[spec.centers_lat[i, j, p] for (i, j, p) in cells]
        lons_center = Float32[spec.centers_lon[i, j, p] for (i, j, p) in cells]
        signal = Float32[10p for (_, _, p) in cells]

        lat_bnds = zeros(Float32, n, 4)
        lon_bnds = zeros(Float32, n, 4)
        dx = 0.05f0
        for k in 1:n
            lat_bnds[k, :] = [lats_center[k] - dx, lats_center[k] - dx,
                              lats_center[k] + dx, lats_center[k] + dx]
            lon_bnds[k, :] = [lons_center[k] - dx, lons_center[k] + dx,
                              lons_center[k] + dx, lons_center[k] - dx]
        end

        input_file = joinpath(test_dir, "TROPO_SIF_2020-01-01_cs.nc")
        ds = Dataset(input_file, "c")
        defDim(ds, "x", n)
        defDim(ds, "nv", 4)
        defVar(ds, "lat", Float32, ("x",))[:] = lats_center
        defVar(ds, "lon", Float32, ("x",))[:] = lons_center
        defVar(ds, "lat_bnds", Float32, ("x", "nv"))[:, :] = lat_bnds
        defVar(ds, "lon_bnds", Float32, ("x", "nv"))[:, :] = lon_bnds
        defVar(ds, "sif", Float32, ("x",))[:] = signal
        defVar(ds, "sza", Float32, ("x",))[:] = fill(30.0f0, n)
        close(ds)

        config_file = joinpath(test_dir, "cs_config.toml")
        write(config_file, """
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
        """)

        output_file = joinpath(test_dir, "cs_output.nc")
        config = load_config(config_file)
        time_spec = TimeSpec(DateTime("2020-01-01"), DateTime("2020-01-01"), Dates.Day(1))
        grid(config, spec, time_spec, SubpixelGridding(); outfile=output_file)

        ds_out = Dataset(output_file)
        @test size(ds_out["lons"]) == (Nc, Nc, 6)
        @test size(ds_out["corner_lons"]) == (Nc + 1, Nc + 1, 6)
        weights = ds_out["n"][:, :, :, 1]
        sif_grid = ds_out["sif"][:, :, :, 1]
        @test size(weights) == (Nc, Nc, 6, )

        # Conservation: each footprint deposited total weight 1.
        @test sum(weights) ≈ Float32(n) rtol = 1f-4
        # Each target cell holds exactly its footprint's signal.
        for (k, (i, j, p)) in enumerate(cells)
            @test weights[i, j, p] ≈ 1.0f0 atol = 1f-4
            @test sif_grid[i, j, p] ≈ signal[k] rtol = 1f-4
        end
        # Only the n target cells hold data; untouched cells are fill (read as `missing`).
        @test count(!ismissing, sif_grid) == n
        close(ds_out)

        rm(test_dir; recursive=true)
    end

    @testset "Footprint straddling a panel edge deposits on two panels" begin
        test_dir = mktempdir()
        spec = CubedSphereGridSpec(Nc=12, T=Float64)
        Nc = spec.Nc

        # Find an equatorial longitude where the panel changes within a few
        # degrees, then center a wide footprint across that seam.
        panel_at(lon) = lonlat_to_panel_xy(spec.definition, Nc, lon, 0.0, Float64)[1]
        seam = nothing
        lon = 0.0
        while lon < 360.0
            if panel_at(lon) != panel_at(lon + 0.5)
                seam = lon + 0.25
                break
            end
            lon += 0.5
        end
        @test seam !== nothing

        half_w = 3.0  # ±3° in lon spans the seam
        half_h = 1.5
        lat_bnds = reshape(Float64[-half_h, -half_h, half_h, half_h], 1, 4)
        lon_bnds = reshape(Float64[seam - half_w, seam + half_w,
                                   seam + half_w, seam - half_w], 1, 4)

        input_file = joinpath(test_dir, "SEAM_2020-01-01.nc")
        ds = Dataset(input_file, "c")
        defDim(ds, "x", 1)
        defDim(ds, "nv", 4)
        defVar(ds, "lat_bnds", Float64, ("x", "nv"))[:, :] = lat_bnds
        defVar(ds, "lon_bnds", Float64, ("x", "nv"))[:, :] = lon_bnds
        defVar(ds, "sif", Float64, ("x",))[:] = Float64[10.0]
        close(ds)

        config_file = joinpath(test_dir, "seam_config.toml")
        write(config_file, """
        filePattern = "SEAM_YYYY-MM-DD.nc"
        folder = "$(replace(test_dir, "\\" => "/"))/"

        [basic]
        lat_bnd = "lat_bnds"
        lon_bnd = "lon_bnds"

        [grid]
        sif = "sif"
        """)

        output_file = joinpath(test_dir, "seam_output.nc")
        config = load_config(config_file)
        time_spec = TimeSpec(DateTime("2020-01-01"), DateTime("2020-01-01"), Dates.Day(1))
        # Force oversampling so the seam-crossing subdivision is exercised.
        grid(config, spec, time_spec, SubpixelGridding(n_oversample=10); outfile=output_file)

        ds_out = Dataset(output_file)
        weights = ds_out["n"][:, :, :, 1]
        sif_grid = ds_out["sif"][:, :, :, 1]
        # Conservation across the seam.
        @test sum(weights) ≈ 1.0 rtol = 1e-3
        # Weight landed on at least two distinct panels.
        panels_hit = count(p -> any(weights[:, :, p] .> 0), 1:6)
        @test panels_hit >= 2
        # Every populated cell carries the single footprint's value.
        for idx in CartesianIndices(sif_grid)
            if weights[idx] > 1e-10 && sif_grid[idx] > -999
                @test sif_grid[idx] ≈ 10.0 rtol = 1e-6
            end
        end
        close(ds_out)

        rm(test_dir; recursive=true)
    end

    @testset "Unsupported combinations error clearly" begin
        spec = CubedSphereGridSpec(Nc=6)
        config = DataSourceConfig(Dict("lat_bnd" => "lat_bnds", "lon_bnd" => "lon_bnds"),
                                  OrderedDict("sif" => "sif"), FilterRule[],
                                  "X_YYYY-MM-DD.nc", "/tmp/")
        time_spec = TimeSpec(DateTime("2020-01-01"), DateTime("2020-01-01"), Dates.Day(1))
        # GPU/KA backend not supported on the cube yet.
        @test_throws ErrorException grid(config, spec, time_spec, SubpixelGridding();
                                         backend=CPU())
        # Circular footprints not supported on the cube yet.
        @test_throws ErrorException grid(config, spec, time_spec,
                                         CircularFootprintGridding())
    end
end
