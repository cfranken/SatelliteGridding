function _circular_accumulation_allocations(n_footprints)
    means = zeros(Float32, 2, 2, 1)
    weights = zeros(Float32, 2, 2)
    centers = fill(1.5f0, n_footprints)
    corners = repeat(Float32[1.2 1.2 1.8 1.8], n_footprints, 1)
    values = fill(1f0, n_footprints, 1)
    accumulate!() = accumulate_circular_footprint!(means, nothing, weights, false,
        centers, centers, corners, corners, values, n_footprints, 1, 20)

    accumulate!()
    @test weights[1, 1] == n_footprints
    @allocated accumulate!()
end

@testset "Averaging" begin

    @testset "Circular allocations are independent of footprint count" begin
        @test _circular_accumulation_allocations(1000) == _circular_accumulation_allocations(10)
    end

    @testset "Mean-only accumulation needs no variance workspace" begin
        for T in (Float32, Float64)
            means = zeros(T, 3, 3, 1)
            weights = zeros(T, 3, 3)
            values = reshape(T[10, 30], 2, 1)
            corners = fill(T(1.5), 2, 4)
            n = 4
            accumulate_footprint!(means, nothing, weights, false, corners, corners,
                values, 2, 1, n, zeros(T, n, n, 2), zeros(Int32, n^2), zeros(Int32, n^2),
                zeros(T, n), zeros(T, n), zeros(T, n), zeros(T, n))
            @test means[1, 1, 1] ≈ 20
            @test weights[1, 1] == 2
        end
    end

    @testset "Welford single-cell — fast path (all corners in one cell)" begin
        # Setup: 5×5 grid, 1 variable
        grid_data = zeros(Float32, 5, 5, 1)
        grid_std = zeros(Float32, 5, 5, 1)
        grid_weights = zeros(Float32, 5, 5)

        n_over = 10
        points = zeros(Float32, n_over, n_over, 2)
        ix = zeros(Int32, n_over^2)
        iy = zeros(Int32, n_over^2)
        lats0 = zeros(n_over); lons0 = zeros(n_over)
        lats1 = zeros(n_over); lons1 = zeros(n_over)

        # 100 observations all with corners inside cell (3, 3)
        # Corners at fractional indices: all within [3.1, 3.9]
        n_pixels = 100
        values = randn(Float32, n_pixels, 1) .+ 10.0f0  # mean ≈ 10
        lat_idx = fill(Float32(3.5), n_pixels, 4)  # all corners same cell
        lon_idx = fill(Float32(3.5), n_pixels, 4)

        accumulate_footprint!(grid_data, grid_std, grid_weights, true,
                              lat_idx, lon_idx, values, n_pixels, 1, n_over,
                              points, ix, iy, lats0, lons0, lats1, lons1)

        # Check Welford mean matches Statistics.mean
        @test grid_data[3, 3, 1] ≈ mean(values[:, 1]) rtol=1e-5
        @test grid_weights[3, 3] ≈ Float32(n_pixels)

        # Check Welford std matches Statistics.std (population std)
        welford_std = sqrt(grid_std[3, 3, 1] / grid_weights[3, 3])
        expected_std = std(values[:, 1]; corrected=false)
        @test welford_std ≈ expected_std rtol=1e-4

        # All other cells should be empty
        @test sum(grid_weights) ≈ grid_weights[3, 3]
    end

    @testset "Oversampling — footprint spanning 4 cells equally" begin
        # Setup: 10×10 grid, 1 variable
        grid_data = zeros(Float32, 10, 10, 1)
        grid_std = zeros(Float32, 10, 10, 1)
        grid_weights = zeros(Float32, 10, 10)

        n_over = 10
        points = zeros(Float32, n_over, n_over, 2)
        ix = zeros(Int32, n_over^2)
        iy = zeros(Int32, n_over^2)
        lats0 = zeros(n_over); lons0 = zeros(n_over)
        lats1 = zeros(n_over); lons1 = zeros(n_over)

        # Single footprint centered at grid cell intersection (5.0, 5.0)
        # spanning cells (4,4), (4,5), (5,4), (5,5)
        # Corner indices: approximately 4.0 to 6.0 in both dimensions
        V = 42.0f0
        n_pixels = 1
        values = fill(V, 1, 1)

        # Corners define a 2-unit-wide footprint centered at (5,5)
        lat_idx = Float32[4.0 4.0 6.0 6.0]
        lon_idx = Float32[4.0 6.0 6.0 4.0]

        accumulate_footprint!(grid_data, grid_std, grid_weights, false,
                              lat_idx, lon_idx, values, n_pixels, 1, n_over,
                              points, ix, iy, lats0, lons0, lats1, lons1)

        # Each of the 4 cells should get roughly 25 sub-pixels (= weight 0.25)
        # With n=10 and a 2-unit footprint centered at integer boundary,
        # sub-pixels distribute across cells [4,4], [4,5], [5,4], [5,5]
        total_weight = sum(grid_weights)
        @test total_weight ≈ 1.0f0 rtol=1e-5  # Total weight = 1 (one footprint)

        # Each quadrant should get ≈ 0.25 weight
        for ci in 4:5, cj in 4:5
            @test grid_weights[ci, cj] > 0.1  # Each cell gets substantial weight
        end

        # All 4 cells should have value V (since same value everywhere)
        for ci in 4:5, cj in 4:5
            if grid_weights[ci, cj] > 0
                @test grid_data[ci, cj, 1] ≈ V rtol=1e-4
            end
        end
    end

    @testset "Oversampling — footprint spanning 2 cells" begin
        grid_data = zeros(Float32, 10, 10, 1)
        grid_std = zeros(Float32, 10, 10, 1)
        grid_weights = zeros(Float32, 10, 10)

        n_over = 10
        points = zeros(Float32, n_over, n_over, 2)
        ix = zeros(Int32, n_over^2)
        iy = zeros(Int32, n_over^2)
        lats0 = zeros(n_over); lons0 = zeros(n_over)
        lats1 = zeros(n_over); lons1 = zeros(n_over)

        # Footprint: 1 unit tall, 2 units wide, straddling lon boundary at 5.0
        # Spans cells (5,4) and (5,5) in lon
        V = 7.0f0
        values = fill(V, 1, 1)
        lat_idx = Float32[5.1 5.1 5.9 5.9]  # all in lat cell 5
        lon_idx = Float32[4.0 6.0 6.0 4.0]  # spans lon cells 4-5

        accumulate_footprint!(grid_data, grid_std, grid_weights, false,
                              lat_idx, lon_idx, values, 1, 1, n_over,
                              points, ix, iy, lats0, lons0, lats1, lons1)

        # Total weight should be 1
        @test sum(grid_weights) ≈ 1.0f0 rtol=1e-5

        # Weight should be split approximately 50/50 across the two lon cells
        w4 = grid_weights[4, 5] + grid_weights[5, 5]  # two possible lon cells for lat=5
        @test w4 > 0.8  # Most weight in these two cells

        # Value should be V everywhere it has weight
        for i in 1:10, j in 1:10
            if grid_weights[i, j] > 0.01
                @test grid_data[i, j, 1] ≈ V rtol=1e-3
            end
        end
    end

    @testset "Multiple footprints — Welford mean and std" begin
        grid_data = zeros(Float32, 5, 5, 1)
        grid_std = zeros(Float32, 5, 5, 1)
        grid_weights = zeros(Float32, 5, 5)

        n_over = 10
        points = zeros(Float32, n_over, n_over, 2)
        ix = zeros(Int32, n_over^2)
        iy = zeros(Int32, n_over^2)
        lats0 = zeros(n_over); lons0 = zeros(n_over)
        lats1 = zeros(n_over); lons1 = zeros(n_over)

        # Two footprints, both fully inside cell (3,3) (fast path)
        V1, V2 = 10.0f0, 20.0f0

        for V in [V1, V2]
            values = fill(V, 1, 1)
            lat_idx = fill(Float32(3.5), 1, 4)
            lon_idx = fill(Float32(3.5), 1, 4)
            accumulate_footprint!(grid_data, grid_std, grid_weights, true,
                                  lat_idx, lon_idx, values, 1, 1, n_over,
                                  points, ix, iy, lats0, lons0, lats1, lons1)
        end

        @test grid_data[3, 3, 1] ≈ (V1 + V2) / 2 rtol=1e-5
        @test grid_weights[3, 3] ≈ 2.0f0

        # Welford std = population std of [10, 20] = 5.0
        welford_std = sqrt(grid_std[3, 3, 1] / grid_weights[3, 3])
        @test welford_std ≈ std([V1, V2]; corrected=false) rtol=1e-4
    end

    @testset "No cross-contamination between cells" begin
        grid_data = zeros(Float32, 10, 10, 1)
        grid_std = zeros(Float32, 10, 10, 1)
        grid_weights = zeros(Float32, 10, 10)

        n_over = 10
        points = zeros(Float32, n_over, n_over, 2)
        ix = zeros(Int32, n_over^2)
        iy = zeros(Int32, n_over^2)
        lats0 = zeros(n_over); lons0 = zeros(n_over)
        lats1 = zeros(n_over); lons1 = zeros(n_over)

        # Footprint A in cell (2,2), value=100
        values_a = fill(100.0f0, 1, 1)
        lat_a = fill(Float32(2.5), 1, 4)
        lon_a = fill(Float32(2.5), 1, 4)
        accumulate_footprint!(grid_data, grid_std, grid_weights, false,
                              lat_a, lon_a, values_a, 1, 1, n_over,
                              points, ix, iy, lats0, lons0, lats1, lons1)

        # Footprint B in cell (8,8), value=200
        values_b = fill(200.0f0, 1, 1)
        lat_b = fill(Float32(8.5), 1, 4)
        lon_b = fill(Float32(8.5), 1, 4)
        accumulate_footprint!(grid_data, grid_std, grid_weights, false,
                              lat_b, lon_b, values_b, 1, 1, n_over,
                              points, ix, iy, lats0, lons0, lats1, lons1)

        @test grid_data[2, 2, 1] ≈ 100.0f0
        @test grid_data[8, 8, 1] ≈ 200.0f0
        @test grid_weights[2, 2] ≈ 1.0f0
        @test grid_weights[8, 8] ≈ 1.0f0
        # Everything else should be zero
        @test sum(grid_weights) ≈ 2.0f0
    end

    @testset "Welford numerical stability" begin
        # Values near 1e8 with small perturbations — naive averaging would lose precision
        grid_data = zeros(Float64, 1, 1, 1)
        grid_std = zeros(Float64, 1, 1, 1)
        grid_weights = zeros(Float64, 1, 1)

        n_over = 10
        points = zeros(Float64, n_over, n_over, 2)
        ix = zeros(Int32, n_over^2)
        iy = zeros(Int32, n_over^2)
        lats0 = zeros(n_over); lons0 = zeros(n_over)
        lats1 = zeros(n_over); lons1 = zeros(n_over)

        base_val = 1e8
        perturbations = randn(1000) .* 0.01
        all_values = base_val .+ perturbations

        for v in all_values
            values = fill(v, 1, 1)
            lat_idx = fill(1.5, 1, 4)
            lon_idx = fill(1.5, 1, 4)
            accumulate_footprint!(grid_data, grid_std, grid_weights, true,
                                  lat_idx, lon_idx, values, 1, 1, n_over,
                                  points, ix, iy, lats0, lons0, lats1, lons1)
        end

        @test grid_data[1, 1, 1] ≈ mean(all_values) rtol=1e-10
        welford_std = sqrt(grid_std[1, 1, 1] / grid_weights[1, 1])
        @test welford_std ≈ std(all_values; corrected=false) rtol=1e-6
    end

    @testset "accumulate_center!" begin
        grid_data = zeros(Float32, 5, 5, 2)
        grid_weights = zeros(Float32, 5, 5)

        lat_idx = Int32[2, 2, 4]
        lon_idx = Int32[3, 3, 1]
        values = Float32[10 20; 30 40; 50 60]

        accumulate_center!(grid_data, grid_weights, lat_idx, lon_idx, values, 3, 2)

        @test grid_weights[3, 2] ≈ 2.0f0  # two pixels in cell (3,2)
        @test grid_weights[1, 4] ≈ 1.0f0  # one pixel in cell (1,4)
        @test grid_data[3, 2, 1] ≈ 40.0f0  # sum of 10+30
        @test grid_data[3, 2, 2] ≈ 60.0f0  # sum of 20+40
        @test grid_data[1, 4, 1] ≈ 50.0f0
    end

    @testset "finalize_std!" begin
        grid_std = Float32[4.0;;; 9.0]  # M2 values, shape (1,1,2)
        grid_weights = Float32[4.0;;]   # shape (1,1)

        SatelliteGridding.finalize_std!(grid_std, grid_weights)

        @test grid_std[1, 1, 1] ≈ 1.0f0  # sqrt(4/4)
        @test grid_std[1, 1, 2] ≈ sqrt(9.0f0 / 4.0f0)
    end

    @testset "accumulate_footprint! bounds guard — prevfloat(90f0) through grid_l2" begin
        spec = GridSpec(lat_min=-90f0, lat_max=90f0,
                        lon_min=-180f0, lon_max=180f0,
                        dlat=2f0, dlon=2f0)
        edge_lat = prevfloat(90f0)  # passes apply_filters (< 90f0) but maps to row 91.0f0 in Float32
        mktempdir() do dir
            infile = joinpath(dir, "edge_seq_2020-01-01.nc")
            Dataset(infile, "c") do ds
                defDim(ds, "x", 3)
                defDim(ds, "nv", 4)
                # 1: fast-path all corners at prevfloat(90f0) -> row 91 (OOB, dropped)
                # 2: oversample-path spanning [87.5f0, prevfloat(90f0)] -> partial in-grid
                # 3: interior footprint in cell (91, 46)
                defVar(ds, "lat", Float32, ("x",))[:] = Float32[edge_lat, 88.75f0, 0.5f0]
                defVar(ds, "lon", Float32, ("x",))[:] = Float32[0.5f0, 4.5f0, 0.5f0]
                defVar(ds, "lat_bnds", Float32, ("x", "nv"))[:, :] = Float32[
                    edge_lat edge_lat edge_lat edge_lat;
                    87.5f0   87.5f0   edge_lat edge_lat;
                    0.2f0    0.2f0    0.8f0    0.8f0
                ]
                defVar(ds, "lon_bnds", Float32, ("x", "nv"))[:, :] = Float32[
                    0.2f0 0.8f0 0.8f0 0.2f0;
                    4.2f0 4.8f0 4.8f0 4.2f0;
                    0.2f0 0.8f0 0.8f0 0.2f0
                ]
                defVar(ds, "val", Float32, ("x",))[:] = Float32[99.0f0, 20.0f0, 5.0f0]
            end
            cfg = DataSourceConfig(
                Dict("lat" => "lat", "lon" => "lon",
                     "lat_bnd" => "lat_bnds", "lon_bnd" => "lon_bnds"),
                OrderedDict("val" => "val"), FilterRule[], "edge_seq_YYYY-MM-DD.nc", dir)
            ts = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 1), Dates.Day(1))
            outfile = joinpath(dir, "out_seq.nc")

            grid_l2(cfg, spec, ts; footprint_method=SubpixelGridding(4),
                    outfile=outfile, backend=nothing, progress=false)

            w, v = Dataset(outfile) do ds
                Array(ds["n"][:, :, 1]), Array(ds["val"][:, :, 1])
            end
            # Fast-path prevfloat(90f0) footprint was dropped without OOB write
            @test w[91, 90] == 0.0f0
            # Interior footprint accumulated normally
            @test w[91, 46] == 1.0f0
            @test v[91, 46] ≈ 5.0f0
            # Oversample footprint landed cleanly in column 93 (lon 4..6°), rows 89..90
            @test w[93, 89] + w[93, 90] > 0.0f0
            @test w[93, 89] + w[93, 90] <= 1.0f0
            @test v[93, 90] ≈ 20.0f0
        end
    end
end
