using SHA, Random

function _cubed_sphere_footprints(::Type{T}) where {T}
    # Repeated, overlapping footprints across a panel edge and the dateline,
    # plus a near-pole footprint and one wholly inside a cell.
    lat = [-2.2 -2.1 1.9 2.6; -1.4 -1.3 2.0 2.2; 3.0 3.1 7.6 7.0;
           4.0 4.0 8.0 8.0; 86.0 87.0 89.0 89.2; -0.02 -0.02 0.02 0.02]
    lon = [32.1 49.3 47.7 34.0; 33.2 50.0 46.8 37.0; 176.1 -176.2 -177.0 177.4;
           177.2 -178.0 -177.3 177.4; 127.0 -137.0 -101.0 113.0;
           -10.02 -9.98 -9.98 -10.02]
    values = T[1.1 9.2; -4.3 6.4; 13.5 -2.6; 31.7 5.8; 0.9 4.1; 8.2 1.3]
    return lat, lon, values
end

function _cubed_sphere_accumulator_bits(spec::CubedSphereGridSpec{T}) where {T}
    lat, lon, values = _cubed_sphere_footprints(T)
    n = 7
    shape = accumulator_shape(spec)
    data, m2, weights = zeros(T, shape..., 2), zeros(T, shape..., 2), zeros(T, shape)
    observations = (lat_bounds=lat, lon_bounds=lon, values=values)
    SatelliteGridding._accumulate_l2_spherical!(
        data, m2, weights, true, spec, observations, n)
    # Serialize the raw bits in a fixed byte order before hashing.
    bits_type = T == Float32 ? UInt32 : UInt64
    bits = reinterpret(bits_type, vcat(vec(data), vec(m2), vec(weights)))
    bytes2hex(sha256(reinterpret(UInt8, htol.(bits))))
end

function _cube_reference_cell(def, Nc, lon, lat)
    panel, s, t = lonlat_to_panel_xy(def, Nc, lon, lat, Float64)
    return floor(Int, s) + Nc * (floor(Int, t) - 1 + Nc * (panel - 1))
end

function _cube_footprint_weights(spec, lat, lon, n)
    shape = accumulator_shape(spec)
    weights = zeros(shape)
    lat_bounds = reshape(collect(lat), 1, 4)
    lon_bounds = reshape(collect(lon), 1, 4)
    SatelliteGridding._order_l2_corners!(spec, lat_bounds, lon_bounds)
    observations = (; lat_bounds, lon_bounds, values=ones(1, 1))
    SatelliteGridding._accumulate_l2_spherical!(
        zeros(shape..., 1), nothing, weights, false, spec, observations, n)
    return weights
end

# Compare occupied cells using the four bilinear basis functions and the
# public inverse projection, independently of production's tangent sort.
function _cube_dense_reference(spec, vertices, n)
    cells = Set{Int}()
    for j in 1:n, i in 1:n
        u, v = (i - 0.5) / n, (j - 0.5) / n
        basis = ((1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v)
        x, y, z = ntuple(axis -> sum(basis[k] * vertices[k][axis] for k in 1:4), 3)
        push!(cells, _cube_reference_cell(spec.definition, spec.Nc,
                                          atand(y, x), atand(z, hypot(x, y))))
    end
    return cells
end

function _test_cube_ka(backend)
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        lat, lon, values = _cubed_sphere_footprints(Float64)
        for day in 1:3
            Dataset(joinpath(dir, "cube_2020-01-0$day.nc"), "c") do ds
                defDim(ds, "obs", size(values, 1))
                defDim(ds, "vertex", 4)
                defVar(ds, "lat", Float64, ("obs",))[:] = lat[:, 1]
                defVar(ds, "lon", Float64, ("obs",))[:] = lon[:, 1]
                defVar(ds, "lat_bounds", Float64, ("obs", "vertex"))[:, :] = lat
                defVar(ds, "lon_bounds", Float64, ("obs", "vertex"))[:, :] = lon
                defVar(ds, "signal", Float64, ("obs",))[:] = values[:, 1] .+ day
            end
        end
        config = DataSourceConfig(
            Dict("lat_bnd" => "lat_bounds", "lon_bnd" => "lon_bounds",
                 "lat" => "lat", "lon" => "lon"),
            OrderedDict("signal" => "signal"), FilterRule[], "cube_YYYY-MM-DD.nc", dir)
        config.options["output"] = Dict("precision" => "float64",
            "significant_digits_mean" => "none", "significant_digits_std" => "none")
        start, stop = DateTime(2020, 1, 1), DateTime(2020, 1, 2)
        tiled = TimeSpec(start, stop, Day(1))
        rolling = TimeSpec(start, stop, Day(1); sample_step=Day(1), window_halfwidth_days=1)
        seq_path, ka_path = joinpath(dir, "seq.nc"), joinpath(dir, "ka.nc")
        for def in (GMAOCubedSphereDefinition(), EquiangularCubedSphereDefinition()),
                (times, n, compute_std, expected) in (
                    (tiled, 7, true, [6, 6]), (rolling, nothing, true, [12, 18]),
                    (tiled, nothing, false, [6, 6]))
            spec = CubedSphereGridSpec(Nc=9, T=Float64, definition=def)
            opts = (; n_oversample=n, compute_std, progress=false, verbose=false)
            grid_l2(config, spec, times; outfile=seq_path, opts...)
            for (geometry_precision, accumulator_precision) in ((nothing, nothing), (Float32, Float32))
                grid_l2(config, spec, times; backend, geometry_precision, accumulator_precision,
                        outfile=ka_path, opts...)
                Dataset(seq_path) do seq
                    Dataset(ka_path) do ka
                        w_seq, w_ka = seq["n"][:, :, :, :], ka["n"][:, :, :, :]
                        hits = w_seq .> 0
                        @test hits == (w_ka .> 0)
                        tol = geometry_precision === Float32 ? 2e-6 : 1e-12
                        if backend isa CPU && geometry_precision === nothing
                            @test w_seq == w_ka
                        else
                            @test w_seq ≈ w_ka rtol=tol
                        end
                        @test vec(sum(w_ka; dims=(1, 2, 3))) ≈ expected rtol=tol
                        @test all(isapprox.(seq["signal"][:, :, :, :][hits],
                                          ka["signal"][:, :, :, :][hits]; rtol=tol, atol=tol))
                        @test haskey(ka, "signal_std") == compute_std
                        if compute_std
                            @test all(isapprox.(seq["signal_std"][:, :, :, :][hits],
                                              ka["signal_std"][:, :, :, :][hits]; rtol=tol, atol=tol))
                        end
                    end
                end
            end
        end
    end
end

@testset "Cubed-sphere grid" begin

    @testset "Unit-vector sequential accumulation bits" begin
        # SHA-256 of the mean, M2, and weight bits. The 7×7 samples use non-binary weights.
        expected = (
            ("2875fd8d716eee69725037219282d3ff66f081e1b2135f049298d49bed0c360e",
             "80db235dfa84f380a5a29cca9a5525263b6012414d694b9bc71c3836ac5df644"),
            ("139005b60a9c430b0abee9f3e2f604ba3cc7c4e164ca9c565adfa79c9f5dee5d",
             "07808d5213b2329f3c9e703897b05c3f29968cd4ff43d15f712a8ee64ea16e23"),
        )
        definitions = (GMAOCubedSphereDefinition(), EquiangularCubedSphereDefinition())
        for (d, def) in enumerate(definitions), (t, T) in enumerate((Float32, Float64))
            spec = CubedSphereGridSpec(Nc=9, definition=def, T=T)
            @test_bitwise _cubed_sphere_accumulator_bits(spec) == expected[d][t]
        end
    end

    @testset "Device locator precision" begin
        S = SatelliteGridding
        rng = MersenneTwister(782)
        lons, lats = 360 .* rand(rng, 500_000), asind.(2 .* rand(rng, 500_000) .- 1)
        for def in (GMAOCubedSphereDefinition(), EquiangularCubedSphereDefinition())
            coarse = CubedSphereGridSpec(Nc=48, definition=def, T=Float64)
            fine = CubedSphereGridSpec(Nc=720, definition=def, T=Float32)
            signed_corners = [S._corner_xyz(def, coarse.Nc, i, j, panel, Float64)
                              for i in 1:49, j in 1:49, panel in 1:6]
            signed_lon = map(v -> atand(v[2], v[1]), signed_corners)
            signed_lat = map(v -> asind(v[3] / sqrt(sum(abs2, v))), signed_corners)
            for (spec, lon, lat) in ((coarse, coarse.centers_lon, coarse.centers_lat),
                    (coarse, coarse.corners_lon, coarse.corners_lat),
                    (coarse, signed_lon, signed_lat),
                    (fine, lons, lats), (fine, lons .- 180, lats))
                f64, f32 = S._device_locator(spec), S._device_locator(spec, Float32)
                @test isbits(f64) && isbits(f32)
                mismatches64 = mismatches32 = 0
                for (x, y) in zip(lon, lat)
                    id = _cube_reference_cell(def, spec.Nc, x, y)
                    point = S._sphere_unit_vector(spec, y, x)
                    mismatches64 += S.locate(f64, point...) != id
                    mismatches32 += S.locate(f32, Float32.(point)...) != id
                end
                @test mismatches64 == 0
                if spec === fine
                    @test mismatches32 / length(lon) < 1e-4
                elseif lon === coarse.centers_lon
                    @test mismatches32 == 0
                end
            end

            # Bracket the equatorial panel edge on both sides, approximately
            # two metres or more from the boundary.
            edge = S.longitude_offset_deg(def) + 45
            f32 = S._device_locator(fine, Float32)
            @test all(S.locate(f32, Float32.(S._sphere_unit_vector(fine, lat, edge + dx))...) ==
                      _cube_reference_cell(def, fine.Nc, edge + dx, lat)
                      for dx in (-2e-4, -2e-5, 2e-5, 2e-4), lat in (-27.13, 0.19, 23.47))
        end
    end

    @testset "Signed panel-edge longitudes in Float64" begin
        S = SatelliteGridding
        for (def, lon, cell) in ((GMAOCubedSphereDefinition(), prevfloat(-55.0), 11520),
                                 (EquiangularCubedSphereDefinition(), prevfloat(-45.0), 6960))
            spec = CubedSphereGridSpec(Nc=48, definition=def, T=Float64)
            lat = -34.0
            @test _cube_reference_cell(def, spec.Nc, lon, lat) == cell
            @test S.locate(S._device_locator(spec), S._sphere_unit_vector(spec, lat, lon)...) == cell

            # A collapsed quadrilateral samples the signed edge point itself.
            observations = (lat_bounds=fill(lat, 1, 4), lon_bounds=fill(lon, 1, 4), values=ones(1, 1))
            for backend in (nothing, CPU())
                options = S._spherical_ka_options(spec, backend, nothing, nothing)
                acc = S.new_accumulator(spec, backend, options, 1, false)
                S._accumulate!(acc, spec, observations, SubpixelGridding(1))
                weights = backend === nothing ? acc.weights : acc.moments.weights
                @test findall(!iszero, vec(weights)) == [cell]
                @test weights[cell] == 1
            end
        end
    end

    @testset "Seams, corner order, and enclosed poles" begin
        S = SatelliteGridding
        for def in (GMAOCubedSphereDefinition(), EquiangularCubedSphereDefinition())
            spec = CubedSphereGridSpec(Nc=24, definition=def, T=Float64)
            seam_spec = CubedSphereGridSpec(Nc=48, definition=def, T=Float64)
            lat = (-4.0, -4.0, 4.0, 4.0)
            lon = (177.0, 179.0, 181.0, 179.0)
            wrapped = (177.0, 179.0, -179.0, 179.0)
            perm = [3, 1, 4, 2]
            seam = _cube_footprint_weights(seam_spec, lat, lon, 24)
            @test seam == _cube_footprint_weights(seam_spec, lat, wrapped, 24)
            @test seam == _cube_footprint_weights(seam_spec, lat[perm], wrapped[perm], 24)
            @test sum(seam) ≈ 1

            edge = S.longitude_offset_deg(def) + 45
            lat = (-3.0, -3.0, 3.0, 3.0)
            lon = (edge - 3, edge + 3, edge + 3.5, edge - 2.5)
            panel_edge = _cube_footprint_weights(spec, lat, lon, 8)
            @test panel_edge == _cube_footprint_weights(spec, lat[perm], lon[perm] .- 360, 8)
            @test sum(panel_edge) == 1
            @test vec(sum(reshape(panel_edge, spec.Nc, spec.Nc, 6); dims=(1, 2)))[3:6] == zeros(4)
            @test all(>(0), vec(sum(reshape(panel_edge, spec.Nc, spec.Nc, 6); dims=(1, 2)))[1:2])

            for pole in (-89.0, 89.0)
                lat = (pole, pole, pole, pole)
                lon = S.longitude_offset_deg(def) .+ (-135.0, -45.0, 45.0, 135.0)
                polar = _cube_footprint_weights(spec, lat, lon, 8)
                @test filter(!iszero, vec(polar)) == fill(0.25, 4)
                @test sum(polar) == 1
                @test polar == _cube_footprint_weights(spec, lat[perm], lon[perm], 8)
            end
        end
    end

    @testset "Adaptive sampling covers a wide footprint on C360" begin
        S = SatelliteGridding
        for def in (GMAOCubedSphereDefinition(), EquiangularCubedSphereDefinition())
            spec = CubedSphereGridSpec(Nc=360, definition=def, T=Float64)
            vertices = map(((s, t),) -> S._continuous_panel_xyz(def, spec.Nc, s, t, 1, Float64),
                           ((100.24, 130.27), (106.8, 130.27), (106.8, 132.72), (100.24, 132.72)))
            lon = map(v -> atand(v[2], v[1]), vertices)
            lat = map(v -> atand(v[3], hypot(v[1], v[2])), vertices)
            dense = _cube_dense_reference(spec, vertices, 512)
            adaptive = _cube_footprint_weights(spec, lat, lon, nothing)
            fixed = _cube_footprint_weights(spec, lat, lon, 4)
            @test Set(findall(!iszero, vec(adaptive))) == dense
            @test count(!iszero, fixed) < length(dense)
            @test adaptive == _cube_footprint_weights(spec, lat, lon, 24)
            @test sum(adaptive) ≈ 1
        end
    end

    @testset "Deprecated accumulator validates and delegates" begin
        spec = CubedSphereGridSpec(Nc=9)
        shape = accumulator_shape(spec)
        lat, lon, values = _cubed_sphere_footprints(Float32)
        data, weights = zeros(Float32, shape..., 1), zeros(Float32, shape)
        lat[2, 1] = NaN # Rows beyond n_pixels are ignored.
        @test_deprecated accumulate_footprint_cs!(data, nothing, weights, false, spec,
            lat, lon, values, 1, 1, 7, nothing, nothing, nothing, nothing, nothing)
        @test weights == Float32.(_cube_footprint_weights(spec, lat[1, :], lon[1, :], 7))
        @test all(==(values[1, 1]), view(data, :, :, 1)[weights .> 0])

        lat[1, 1] = NaN
        @test_throws ArgumentError accumulate_footprint_cs!(data, nothing, weights, false, spec,
            lat, lon, values, 1, 1, 7, nothing, nothing, nothing, nothing, nothing)
        lat[1, 1] = lat[1, 2]
        for n in (-1, 0)
            @test_throws ArgumentError accumulate_footprint_cs!(data, nothing, weights, false, spec,
                lat, lon, values, 1, 1, n, nothing, nothing, nothing, nothing, nothing)
        end
        @test sum(weights) ≈ 1
    end

    @testset "KA CPU" _test_cube_ka(CPU())

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
            Base.invokelatest(_test_cube_ka, backend)
        else
            @test_skip available
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
        @test resolve_backend("cpu"; grid_spec=spec) isa CPU
        @test_throws ArgumentError("Metal requires Float32 geometry_precision") resolve_backend(
            "metal"; grid_spec=spec, geometry_precision=Float64)
        # Circular footprints not supported on the cube yet.
        @test_throws ErrorException grid(config, spec, time_spec,
                                         CircularFootprintGridding())
    end
end
