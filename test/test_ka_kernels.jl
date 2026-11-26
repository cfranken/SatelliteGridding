if !@isdefined(_MockNoFloat64Backend)
    struct _MockNoFloat64Backend <: KernelAbstractions.Backend end
    KernelAbstractions.supports_float64(::_MockNoFloat64Backend) = false
end

@testset "KA Kernels" begin
    backend = CPU()

    @testset "sort_corners_ccw_ka! matches scalar version" begin
        # Scrambled corners (like OCO-2: SW, NW, NE, SE)
        lat_mat = Float32[-1 1 1 -1; 0 1 0 -1]
        lon_mat = Float32[-1 -1 1 1; -1 0 1 0]

        # Copy for KA version
        lat_ka = copy(lat_mat)
        lon_ka = copy(lon_mat)

        # Scalar version
        sort_corners_ccw!(lat_mat, lon_mat)
        # KA version
        sort_corners_ccw_ka!(backend, lat_ka, lon_ka)

        @test all(isapprox.(lat_ka, lat_mat))
        @test all(isapprox.(lon_ka, lon_mat))
    end

    @testset "sort_corners_ccw_ka! — many soundings" begin
        n = 1000
        lat_mat = randn(Float32, n, 4)
        lon_mat = randn(Float32, n, 4)

        lat_ka = copy(lat_mat)
        lon_ka = copy(lon_mat)

        sort_corners_ccw!(lat_mat, lon_mat)
        sort_corners_ccw_ka!(backend, lat_ka, lon_ka)

        @test all(isapprox.(lat_ka, lat_mat))
        @test all(isapprox.(lon_ka, lon_mat))
    end

    @testset "compute_footprint_indices_ka! — fast path (single cell)" begin
        # All corners in cell (3, 3)
        lat_corners = Float32[3.2 3.3 3.4 3.1]
        lon_corners = Float32[3.5 3.6 3.7 3.4]
        n = 10

        ix = zeros(Int32, 1, n^2)
        iy = zeros(Int32, 1, n^2)
        skip = zeros(Int32, 1)

        compute_footprint_indices_ka!(backend, ix, iy, skip, lat_corners, lon_corners, n)

        @test skip[1] == Int32(1)  # fast path
        @test iy[1, 1] == Int32(3)
        @test ix[1, 1] == Int32(3)
        # Rest should be sentinel
        @test all(ix[1, 2:end] .== Int32(-1))
    end

    @testset "compute_footprint_indices_ka! — oversample path" begin
        # Footprint spanning cells 4-5 in both dimensions
        lat_corners = Float32[4.0 4.0 6.0 6.0]
        lon_corners = Float32[4.0 6.0 6.0 4.0]
        n = 10

        ix = zeros(Int32, 1, n^2)
        iy = zeros(Int32, 1, n^2)
        skip = zeros(Int32, 1)

        compute_footprint_indices_ka!(backend, ix, iy, skip, lat_corners, lon_corners, n)

        @test skip[1] == Int32(0)  # oversample
        # All indices should be 4 or 5
        @test all(4 .<= iy[1, :] .<= 5)
        @test all(4 .<= ix[1, :] .<= 5)
    end

    @testset "compute_footprint_indices_ka! — skip (too wide)" begin
        # Footprint spanning more than n cells in lon
        lat_corners = Float32[3.0 3.0 4.0 4.0]
        lon_corners = Float32[1.0 20.0 20.0 1.0]  # 19 cells wide, n=10
        n = 10

        ix = zeros(Int32, 1, n^2)
        iy = zeros(Int32, 1, n^2)
        skip = zeros(Int32, 1)

        compute_footprint_indices_ka!(backend, ix, iy, skip, lat_corners, lon_corners, n)

        @test skip[1] == Int32(2)  # skip
        @test all(ix[1, :] .== Int32(-1))
    end

    @testset "compute_circular_footprint_indices_ka! — masked circle" begin
        center_lat = Float32[3.0]
        center_lon = Float32[3.0]
        lat_corners = Float32[2.0 2.0 4.0 4.0]
        lon_corners = Float32[2.0 4.0 4.0 2.0]
        n = 10

        ix = zeros(Int32, 1, n^2)
        iy = zeros(Int32, 1, n^2)
        inside_count = zeros(Int32, 1)
        skip = zeros(Int32, 1)

        compute_circular_footprint_indices_ka!(backend, ix, iy, inside_count,
                                               skip, center_lat, center_lon,
                                               lat_corners, lon_corners, n)

        @test skip[1] == Int32(0)
        @test inside_count[1] == Int32(80)
        @test count(x -> x > 0, ix) == 80
        @test all(2 .<= ix[ix .> 0] .<= 3)
        @test all(2 .<= iy[iy .> 0] .<= 3)
    end

    @testset "accumulate_batch! — single cell fast path" begin
        grid_sum = zeros(Float32, 10, 10, 1)
        grid_weights = zeros(Float32, 10, 10)

        # 50 observations all inside cell (5, 5)
        n_fp = 50
        vals = randn(Float32, n_fp, 1) .+ 10f0
        lat_corners = fill(Float32(5.5), n_fp, 4)
        lon_corners = fill(Float32(5.5), n_fp, 4)

        accumulate_batch!(backend, grid_sum, grid_weights,
                          lat_corners, lon_corners, vals, 10, 1)

        @test grid_weights[5, 5] ≈ Float32(n_fp) rtol=1e-5
        # Sum should equal sum of values
        @test grid_sum[5, 5, 1] ≈ sum(vals) rtol=1e-4
        # After finalize_mean!, should equal mean
        finalize_mean!(grid_sum, grid_weights)
        @test grid_sum[5, 5, 1] ≈ mean(vals) rtol=1e-4
    end

    @testset "accumulate_batch! — multi-cell footprint" begin
        grid_sum = zeros(Float32, 10, 10, 1)
        grid_weights = zeros(Float32, 10, 10)

        V = 42.0f0
        vals = fill(V, 1, 1)
        # Footprint spanning cells 4-5 in both dims
        lat_corners = Float32[4.0 4.0 6.0 6.0]
        lon_corners = Float32[4.0 6.0 6.0 4.0]

        accumulate_batch!(backend, grid_sum, grid_weights,
                          lat_corners, lon_corners, vals, 10, 1)

        # Total weight should be ~1
        @test sum(grid_weights) ≈ 1.0f0 rtol=1e-5

        # After finalize, all cells with weight should have value V
        finalize_mean!(grid_sum, grid_weights)
        for i in 1:10, j in 1:10
            if grid_weights[i, j] > 0.01
                @test grid_sum[i, j, 1] ≈ V rtol=1e-3
            end
        end
    end

    @testset "accumulate_circular_batch! — symmetric circle" begin
        grid_sum = zeros(Float32, 6, 6, 1)
        grid_weights = zeros(Float32, 6, 6)
        center_lat = Float32[3.0]
        center_lon = Float32[3.0]
        lat_corners = Float32[2.0 2.0 4.0 4.0]
        lon_corners = Float32[2.0 4.0 4.0 2.0]
        vals = Float32[42.0;;]

        accumulate_circular_batch!(backend, grid_sum, grid_weights,
                                   center_lat, center_lon,
                                   lat_corners, lon_corners, vals, 100, 1)

        @test sum(grid_weights) ≈ 1.0f0 rtol=2f-5
        @test grid_weights[2, 2] ≈ 0.25f0 atol=0.01f0
        @test grid_weights[3, 2] ≈ 0.25f0 atol=0.01f0
        @test grid_weights[2, 3] ≈ 0.25f0 atol=0.01f0
        @test grid_weights[3, 3] ≈ 0.25f0 atol=0.01f0

        finalize_mean!(grid_sum, grid_weights)
        for idx in CartesianIndices(grid_weights)
            if grid_weights[idx] > 0
                @test grid_sum[idx[1], idx[2], 1] ≈ 42.0f0
            end
        end
    end

    @testset "accumulate_circular_batch! matches sequential circular path" begin
        n_over = 30
        n_fp = 80
        n_vars = 2

        center_lat = rand(Float32, n_fp) .* 4 .+ 4
        center_lon = rand(Float32, n_fp) .* 4 .+ 4
        radius_lat = rand(Float32, n_fp) .* 0.6f0 .+ 0.3f0
        radius_lon = rand(Float32, n_fp) .* 0.6f0 .+ 0.3f0
        lat_corners = hcat(center_lat .- radius_lat, center_lat .- radius_lat,
                           center_lat .+ radius_lat, center_lat .+ radius_lat)
        lon_corners = hcat(center_lon .- radius_lon, center_lon .+ radius_lon,
                           center_lon .+ radius_lon, center_lon .- radius_lon)
        vals = randn(Float32, n_fp, n_vars)

        grid_sum = zeros(Float32, 12, 12, n_vars)
        grid_w_ka = zeros(Float32, 12, 12)
        accumulate_circular_batch!(backend, grid_sum, grid_w_ka,
                                   center_lat, center_lon,
                                   lat_corners, lon_corners, vals, n_over, n_vars)
        finalize_mean!(grid_sum, grid_w_ka)

        grid_data = zeros(Float32, 12, 12, n_vars)
        grid_std = zeros(Float32, 12, 12, n_vars)
        grid_w_seq = zeros(Float32, 12, 12)
        accumulate_circular_footprint!(grid_data, grid_std, grid_w_seq, false,
                                       center_lat, center_lon,
                                       lat_corners, lon_corners, vals,
                                       n_fp, n_vars, n_over)

        @test (grid_w_ka .> 0) == (grid_w_seq .> 0)
        @test all(isapprox.(grid_w_ka, grid_w_seq; rtol=1e-4, atol=1e-5))
        mask = grid_w_seq .> 0.01f0
        for z in 1:n_vars
            @test all(isapprox.(grid_sum[:, :, z][mask], grid_data[:, :, z][mask];
                                rtol=1e-3, atol=1e-4))
        end
    end

    @testset "accumulate_batch! matches accumulate_footprint!" begin
        # Compare KA batch pipeline against sequential Welford for mean
        n_over = 10
        n_fp = 200
        n_vars = 3

        # Random footprints well inside the 15×15 grid (indices stay in 1..14)
        lat_base = rand(Float32, n_fp) .* 6 .+ 4  # 4..10
        lon_base = rand(Float32, n_fp) .* 6 .+ 4
        extent = rand(Float32, n_fp) .* 1.0f0 .+ 0.1f0  # 0.1..1.1

        lat_corners = hcat(lat_base .- extent, lat_base .- extent,
                           lat_base .+ extent, lat_base .+ extent)
        lon_corners = hcat(lon_base .- extent, lon_base .+ extent,
                           lon_base .+ extent, lon_base .- extent)
        vals = randn(Float32, n_fp, n_vars)

        # --- KA path: sum-based ---
        grid_sum = zeros(Float32, 15, 15, n_vars)
        grid_w_ka = zeros(Float32, 15, 15)
        lat_ka = copy(lat_corners)
        lon_ka = copy(lon_corners)
        accumulate_batch!(backend, grid_sum, grid_w_ka,
                          lat_ka, lon_ka, vals, n_over, n_vars)
        finalize_mean!(grid_sum, grid_w_ka)

        # --- Sequential path: Welford ---
        grid_data = zeros(Float32, 15, 15, n_vars)
        grid_std = zeros(Float32, 15, 15, n_vars)
        grid_w_seq = zeros(Float32, 15, 15)
        points = zeros(Float32, n_over, n_over, 2)
        ix = zeros(Int32, n_over^2)
        iy = zeros(Int32, n_over^2)
        lats0 = zeros(n_over); lons0 = zeros(n_over)
        lats1 = zeros(n_over); lons1 = zeros(n_over)

        # Sort corners the same way KA does
        lat_seq = copy(lat_corners)
        lon_seq = copy(lon_corners)
        sort_corners_ccw!(lat_seq, lon_seq)

        accumulate_footprint!(grid_data, grid_std, grid_w_seq, false,
                              lat_seq, lon_seq, vals, n_fp, n_vars, n_over,
                              points, ix, iy, lats0, lons0, lats1, lons1)

        # Weights must match bit-for-bit (identical cell assignments and order)
        @test isequal(grid_w_ka, grid_w_seq)

        # Means differ only by Welford vs sum/weight accumulation order
        mask = grid_w_seq .> 0.01f0
        for z in 1:n_vars
            ka_vals = grid_sum[:, :, z][mask]
            seq_vals = grid_data[:, :, z][mask]
            @test all(isapprox.(ka_vals, seq_vals; rtol=1e-3, atol=1e-4))
        end
    end

    @testset "finalize_mean!" begin
        grid_sum = Float32[10.0;;; 20.0]  # shape (1,1,2)
        grid_weights = Float32[5.0;;]     # shape (1,1)

        finalize_mean!(grid_sum, grid_weights)

        @test grid_sum[1, 1, 1] ≈ 2.0f0   # 10/5
        @test grid_sum[1, 1, 2] ≈ 4.0f0   # 20/5
    end

    @testset "Circular KA path — edge/pole bounds, determinism, and sequential agreement" begin
        ka_backend = resolve_backend("cpu")

        @testset "compute_circular_footprint_indices_ka! clips out-of-range indices" begin
            spec = GridSpec(lat_min=-90f0, lat_max=90f0,
                            lon_min=-180f0, lon_max=180f0,
                            dlat=2f0, dlon=2f0)
            n_lon, n_lat = grid_shape(spec)  # (180, 90)

            # Footprints near lat_max, lat_min, lon_max, lon_min, poles, and corners,
            # including asymmetric corners inside (-90,90)×(-180,180) whose circle extends
            # past the grid boundary and corners that round to n_lat+1 / n_lon+1 in Float32.
            center_lat = Float32[
                90.95,   # near lat_max (north pole)
                1.05,    # near lat_min (south pole)
                45.5,    # near lon_max
                45.5,    # near lon_min
                90.95,   # NE corner
                1.05,    # SW corner
                90.999999, # Float32 rounds to 91.0f0 (fast path OOB)
                0.5,     # completely below lat_min (fast path OOB)
                45.5,    # completely past lon_max (zero-radius OOB)
            ]
            center_lon = Float32[
                90.5,
                90.5,
                180.95,
                1.05,
                180.95,
                1.05,
                180.999999, # Float32 rounds to 181.0f0
                90.5,
                181.5,
            ]
            # Asymmetric corner offsets so all 4 corners lie inside [1, n+1) for rows 1..6,
            # yet radius = max|corner - center| pushes the circle outside [1, n+1).
            lat_corners = Float32[
                89.8  89.8  90.99 90.99;
                1.01  1.01  2.2   2.2;
                45.0  45.0  46.0  46.0;
                45.0  45.0  46.0  46.0;
                89.9  89.9  90.99 90.99;
                1.01  1.01  2.1   2.1;
                91.0  91.0  91.0  91.0;
                0.4   0.4   0.6   0.6;
                45.5  45.5  45.5  45.5;
            ]
            lon_corners = Float32[
                90.0   91.0   91.0   90.0;
                90.0   91.0   91.0   90.0;
                179.8  180.99 180.99 179.8;
                1.01   2.2    2.2    1.01;
                179.9  180.99 180.99 179.9;
                1.01   2.1    2.1    1.01;
                181.0  181.0  181.0  181.0;
                90.4   90.6   90.6   90.4;
                181.5  181.5  181.5  181.5;
            ]

            n_fp = length(center_lat)
            for n in (3, 7)
                ix = zeros(Int32, n_fp, n^2)
                iy = zeros(Int32, n_fp, n^2)
                inside_count = zeros(Int32, n_fp)
                skip = zeros(Int32, n_fp)

                compute_circular_footprint_indices_ka!(
                    ka_backend, ix, iy, inside_count, skip,
                    center_lat, center_lon, lat_corners, lon_corners, n;
                    n_lon=n_lon, n_lat=n_lat)

                # No index outside 1..n_lon × 1..n_lat may reach scatter
                @test all(x -> x == Int32(-1) || (Int32(1) <= x <= Int32(n_lon)), ix)
                @test all(y -> y == Int32(-1) || (Int32(1) <= y <= Int32(n_lat)), iy)
                # Fast-path footprints out of bounds must be marked skip (2), not 1
                @test skip[7] == Int32(2)
                @test skip[8] == Int32(2)
                @test skip[9] == Int32(2)
                for fp in 1:n_fp
                    if skip[fp] == Int32(1)
                        @test Int32(1) <= ix[fp, 1] <= Int32(n_lon)
                        @test Int32(1) <= iy[fp, 1] <= Int32(n_lat)
                    end
                end
            end
        end

        @testset "accumulate_circular_batch! edge/pole agreement with sequential and determinism" begin
            n_lon, n_lat = 10, 10
            n_vars = 2
            # Footprints straddling all 4 edges + corners + interior + out-of-grid fast path
            center_lat = Float32[10.8, 1.2, 5.5, 5.5, 10.7, 1.3, 5.5, 11.2, 0.4]
            center_lon = Float32[5.5,  5.5, 10.8, 1.2, 10.7, 1.3, 5.5, 5.5, 5.5]
            r_lat = Float32[0.8, 0.8, 0.6, 0.6, 0.9, 0.9, 0.7, 0.05, 0.0]
            r_lon = Float32[0.6, 0.6, 0.8, 0.8, 0.9, 0.9, 0.7, 0.05, 0.0]
            lat_corners = hcat(center_lat .- r_lat, center_lat .- r_lat,
                               center_lat .+ r_lat, center_lat .+ r_lat)
            lon_corners = hcat(center_lon .- r_lon, center_lon .+ r_lon,
                               center_lon .+ r_lon, center_lon .- r_lon)
            vals = Float32[10 100; 20 200; 30 300; 40 400; 50 500;
                           60 600; 70 700; 80 800; 90 900]
            n_fp = length(center_lat)
            n_over = 11

            # Sequential reference
            grid_data_seq = zeros(Float32, n_lon, n_lat, n_vars)
            grid_w_seq = zeros(Float32, n_lon, n_lat)
            accumulate_circular_footprint!(grid_data_seq, nothing, grid_w_seq, false,
                                           center_lat, center_lon,
                                           lat_corners, lon_corners, vals,
                                           n_fp, n_vars, n_over)

            # Two KA CPU runs must be bit-identical and agree with sequential
            grid_sum_1 = zeros(Float32, n_lon, n_lat, n_vars)
            grid_w_1 = zeros(Float32, n_lon, n_lat)
            accumulate_circular_batch!(ka_backend, grid_sum_1, grid_w_1,
                                       center_lat, center_lon,
                                       lat_corners, lon_corners, vals, n_over, n_vars)

            grid_sum_2 = zeros(Float32, n_lon, n_lat, n_vars)
            grid_w_2 = zeros(Float32, n_lon, n_lat)
            accumulate_circular_batch!(ka_backend, grid_sum_2, grid_w_2,
                                       center_lat, center_lon,
                                       lat_corners, lon_corners, vals, n_over, n_vars)

            @test isequal(grid_w_1, grid_w_2)
            @test isequal(grid_sum_1, grid_sum_2)

            finalize_mean!(grid_sum_1, grid_w_1)
            @test isequal(grid_w_1, grid_w_seq)
            # Partial-weight check: analytic total across 6 clipped + 1 interior footprints is ~4.57
            @test 4.4f0 < sum(grid_w_1) < 4.7f0
            mask = grid_w_seq .> 0.01f0
            for z in 1:n_vars
                @test all(isapprox.(grid_sum_1[:, :, z][mask], grid_data_seq[:, :, z][mask];
                                    rtol=1e-4, atol=1e-4))
            end
        end

        @testset "Shared boundary cell — interior + clipped edge footprint mean (no renormalization)" begin
            n_lon, n_lat = 10, 10
            v_int, v_edge = 10.0f0, 100.0f0
            vals = Float32[v_int; v_edge;;]

            # Circular: interior footprint in cell (5, 10) + clipped edge footprint overlapping (5, 10)
            n_over_c = 11
            center_lat = Float32[10.5, 11.2]
            center_lon = Float32[5.5,  5.5]
            r_lat = Float32[0.2, 0.8]
            r_lon = Float32[0.2, 0.3]
            lat_c = hcat(center_lat .- r_lat, center_lat .- r_lat,
                         center_lat .+ r_lat, center_lat .+ r_lat)
            lon_c = hcat(center_lon .- r_lon, center_lon .+ r_lon,
                         center_lon .+ r_lon, center_lon .- r_lon)

            gs_c = zeros(Float32, n_lon, n_lat, 1)
            gw_c = zeros(Float32, n_lon, n_lat)
            accumulate_circular_batch!(ka_backend, gs_c, gw_c,
                                       center_lat, center_lon, lat_c, lon_c,
                                       vals, n_over_c, 1)
            finalize_mean!(gs_c, gw_c)

            gd_c_seq = zeros(Float32, n_lon, n_lat, 1)
            gw_c_seq = zeros(Float32, n_lon, n_lat)
            accumulate_circular_footprint!(gd_c_seq, nothing, gw_c_seq, false,
                                           center_lat, center_lon, lat_c, lon_c,
                                           vals, 2, 1, n_over_c)

            w_edge_c = gw_c_seq[5, 10] - 1.0f0
            @test 0.05f0 < w_edge_c < 0.8f0
            expected_mean_c = (1.0f0 * v_int + w_edge_c * v_edge) / (1.0f0 + w_edge_c)
            renorm_mean = (v_int + v_edge) / 2.0f0  # 55.0f0 (wrong)
            @test gw_c[5, 10] ≈ gw_c_seq[5, 10] rtol=1e-6
            @test gs_c[5, 10, 1] ≈ expected_mean_c rtol=1e-5
            @test gd_c_seq[5, 10, 1] ≈ expected_mean_c rtol=1e-5
            @test abs(gs_c[5, 10, 1] - renorm_mean) > 10.0f0

            # Quadrilateral: interior footprint (w=1) + half-clipped edge footprint (w=0.5) in cell (5, 10)
            n_over_q = 10
            qlat = Float32[10.2 10.2 10.6 10.6;
                           10.5 10.5 11.5 11.5]
            qlon = Float32[5.2  5.8  5.8  5.2;
                           5.2  5.8  5.8  5.2]

            gs_q = zeros(Float32, n_lon, n_lat, 1)
            gw_q = zeros(Float32, n_lon, n_lat)
            accumulate_batch!(ka_backend, gs_q, gw_q, copy(qlat), copy(qlon), vals, n_over_q, 1)
            finalize_mean!(gs_q, gw_q)

            gd_q_seq = zeros(Float32, n_lon, n_lat, 1)
            gw_q_seq = zeros(Float32, n_lon, n_lat)
            pts = zeros(Float32, n_over_q, n_over_q, 2)
            ixb = zeros(Int32, n_over_q^2); iyb = zeros(Int32, n_over_q^2)
            l0 = zeros(n_over_q); o0 = zeros(n_over_q)
            l1 = zeros(n_over_q); o1 = zeros(n_over_q)
            accumulate_footprint!(gd_q_seq, nothing, gw_q_seq, false,
                                  qlat, qlon, vals, 2, 1, n_over_q,
                                  pts, ixb, iyb, l0, o0, l1, o1)

            expected_mean_q = (1.0f0 * v_int + 0.5f0 * v_edge) / 1.5f0  # 40.0f0
            @test isequal(gw_q, gw_q_seq)
            @test gw_q[5, 10] ≈ 1.5f0 rtol=1e-5
            @test gw_q_seq[5, 10] ≈ 1.5f0 rtol=1e-5
            @test gs_q[5, 10, 1] ≈ expected_mean_q rtol=1e-5
            @test gd_q_seq[5, 10, 1] ≈ expected_mean_q rtol=1e-5
            @test abs(gs_q[5, 10, 1] - renorm_mean) > 10.0f0
        end

        @testset "KA/sequential parity at deterministic edge cases" begin
            # Internal float-precision selection: Float64 on CPU, Float32 when Float64 unsupported
            @test SatelliteGridding._ka_float64_type(ka_backend) === Float64
            @test SatelliteGridding._ka_float64_type(_MockNoFloat64Backend()) === Float32
            @test SatelliteGridding._ka_supported_float(ka_backend, Float32) === Float32
            @test SatelliteGridding._ka_supported_float(ka_backend, Float64) === Float64
            @test SatelliteGridding._ka_supported_float(_MockNoFloat64Backend(), Float32) === Float32
            @test SatelliteGridding._ka_supported_float(_MockNoFloat64Backend(), Float64) === Float32

            # Case 1: footprint at 89.999994°–89.999998° on 2° global grid, n=2 (circular & quad)
            spec = GridSpec(lat_min=-90f0, lat_max=90f0,
                            lon_min=-180f0, lon_max=180f0,
                            dlat=2f0, dlon=2f0)
            mktempdir() do dir
                infile = joinpath(dir, "edge_2020-01-01.nc")
                Dataset(infile, "c") do ds
                    defDim(ds, "x", 1)
                    defDim(ds, "nv", 4)
                    defVar(ds, "lat", Float64, ("x",))[:] = [89.999996]
                    defVar(ds, "lon", Float64, ("x",))[:] = [0.3]
                    defVar(ds, "lat_bnds", Float64, ("x", "nv"))[:, :] =
                        [89.999994 89.999994 89.999998 89.999998]
                    defVar(ds, "lon_bnds", Float64, ("x", "nv"))[:, :] =
                        [0.1 0.5 0.5 0.1]
                    defVar(ds, "val", Float32, ("x",))[:] = Float32[25.0]
                end
                cfg = DataSourceConfig(
                    Dict("lat" => "lat", "lon" => "lon",
                         "lat_bnd" => "lat_bnds", "lon_bnd" => "lon_bnds"),
                    OrderedDict("val" => "val"), FilterRule[], "edge_YYYY-MM-DD.nc", dir)
                ts = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 1), Dates.Day(1))

                for (label, method) in (("circular", CircularFootprintGridding(2)),
                                        ("quad", SubpixelGridding(2)))
                    out_seq = joinpath(dir, "out_seq_$(label).nc")
                    out_ka = joinpath(dir, "out_ka_$(label).nc")
                    grid_l2(cfg, spec, ts; footprint_method=method,
                            outfile=out_seq, backend=nothing, progress=false)
                    grid_l2(cfg, spec, ts; footprint_method=method,
                            outfile=out_ka, backend=ka_backend, progress=false)
                    w_seq, v_seq = Dataset(out_seq) do ds
                        Array(ds["n"][:, :, 1]), Array(ds["val"][:, :, 1])
                    end
                    w_ka, v_ka = Dataset(out_ka) do ds
                        Array(ds["n"][:, :, 1]), Array(ds["val"][:, :, 1])
                    end
                    @test w_seq[91, 90] == 1.0f0
                    @test w_ka[91, 90] == 1.0f0
                    @test isequal(w_ka, w_seq)
                    @test v_ka[91, 90] == 25.0f0
                    @test v_seq[91, 90] == 25.0f0
                end
            end

            # Case 2: Float32 circular kernel rounding case (90.9375f0 / nextfloat(90.8125f0) / 90.99f0, n=2)
            n_lon, n_lat = 180, 90
            clat = Float32[90.9375f0]
            clon = Float32[90.5f0]
            lat_c = Float32[nextfloat(90.8125f0) nextfloat(90.8125f0) 90.99f0 90.99f0]
            lon_c = Float32[90.4f0 90.6f0 90.6f0 90.4f0]
            vals = Float32[17.0;;]
            n = 2

            ix = zeros(Int32, 1, n^2)
            iy = zeros(Int32, 1, n^2)
            ic = zeros(Int32, 1)
            sk = zeros(Int32, 1)
            compute_circular_footprint_indices_ka!(ka_backend, ix, iy, ic, sk,
                                                   clat, clon, lat_c, lon_c, n;
                                                   n_lon=n_lon, n_lat=n_lat)
            @test ic[1] == Int32(4)
            @test all(==(Int32(90)), ix)
            @test all(==(Int32(90)), iy)

            gs_ka = zeros(Float32, n_lon, n_lat, 1)
            gw_ka = zeros(Float32, n_lon, n_lat)
            accumulate_circular_batch!(ka_backend, gs_ka, gw_ka,
                                       clat, clon, lat_c, lon_c, vals, n, 1)
            finalize_mean!(gs_ka, gw_ka)

            gd_seq = zeros(Float32, n_lon, n_lat, 1)
            gw_seq = zeros(Float32, n_lon, n_lat)
            accumulate_circular_footprint!(gd_seq, nothing, gw_seq, false,
                                           clat, clon, lat_c, lon_c, vals, 1, 1, n)

            @test gw_seq[90, 90] == 1.0f0
            @test gw_ka[90, 90] == 1.0f0
            @test isequal(gw_ka, gw_seq)
            @test gs_ka[90, 90, 1] == gd_seq[90, 90, 1]

            # Case 3a: quad sub-sample Float32 rounding before floor — 87.24998474121094°–88.25° × 0.1°–0.5° at n=2
            # Case 3b: pole oversampling drop — 89.999994°–89.999998° × -0.2°–0.2° at n=2
            mktempdir() do dir
                ts = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 1), Dates.Day(1))
                for (CoordT, tag) in ((Float64, "f64"), (Float32, "f32"))
                    f_split = joinpath(dir, "split_$(tag)_2020-01-01.nc")
                    Dataset(f_split, "c") do ds
                        defDim(ds, "x", 1)
                        defDim(ds, "nv", 4)
                        defVar(ds, "lat", CoordT, ("x",))[:] = CoordT[87.75]
                        defVar(ds, "lon", CoordT, ("x",))[:] = CoordT[0.3]
                        defVar(ds, "lat_bnds", CoordT, ("x", "nv"))[:, :] =
                            CoordT[87.24998474121094 87.24998474121094 88.25 88.25]
                        defVar(ds, "lon_bnds", CoordT, ("x", "nv"))[:, :] =
                            CoordT[0.1 0.5 0.5 0.1]
                        defVar(ds, "val", Float32, ("x",))[:] = Float32[25.0]
                    end
                    cfg_split = DataSourceConfig(
                        Dict("lat" => "lat", "lon" => "lon",
                             "lat_bnd" => "lat_bnds", "lon_bnd" => "lon_bnds"),
                        OrderedDict("val" => "val"), FilterRule[],
                        "split_$(tag)_YYYY-MM-DD.nc", dir)
                    out_seq = joinpath(dir, "split_$(tag)_seq.nc")
                    out_ka = joinpath(dir, "split_$(tag)_ka.nc")
                    grid_l2(cfg_split, spec, ts; footprint_method=SubpixelGridding(2),
                            outfile=out_seq, backend=nothing, progress=false)
                    grid_l2(cfg_split, spec, ts; footprint_method=SubpixelGridding(2),
                            outfile=out_ka, backend=ka_backend, progress=false)
                    w_seq, v_seq = Dataset(out_seq) do ds
                        Array(ds["n"][:, :, 1]), Array(ds["val"][:, :, 1])
                    end
                    w_ka, v_ka = Dataset(out_ka) do ds
                        Array(ds["n"][:, :, 1]), Array(ds["val"][:, :, 1])
                    end
                    @test w_seq[91, 89:90] == Float32[0.5f0, 0.5f0]
                    @test w_ka[91, 89:90] == Float32[0.5f0, 0.5f0]
                    @test isequal(w_ka, w_seq)
                    @test v_seq[91, 89:90] == Float32[25.0f0, 25.0f0]
                    @test v_ka[91, 89:90] == Float32[25.0f0, 25.0f0]
                end

                f_drop = joinpath(dir, "pole_drop_2020-01-01.nc")
                Dataset(f_drop, "c") do ds
                    defDim(ds, "x", 1)
                    defDim(ds, "nv", 4)
                    defVar(ds, "lat", Float64, ("x",))[:] = [89.999996]
                    defVar(ds, "lon", Float64, ("x",))[:] = [0.0]
                    defVar(ds, "lat_bnds", Float64, ("x", "nv"))[:, :] =
                        [89.999994 89.999994 89.999998 89.999998]
                    defVar(ds, "lon_bnds", Float64, ("x", "nv"))[:, :] =
                        [-0.2 0.2 0.2 -0.2]
                    defVar(ds, "val", Float32, ("x",))[:] = Float32[25.0]
                end
                cfg_drop = DataSourceConfig(
                    Dict("lat" => "lat", "lon" => "lon",
                         "lat_bnd" => "lat_bnds", "lon_bnd" => "lon_bnds"),
                    OrderedDict("val" => "val"), FilterRule[],
                    "pole_drop_YYYY-MM-DD.nc", dir)
                out_drop_seq = joinpath(dir, "pole_drop_seq.nc")
                out_drop_ka = joinpath(dir, "pole_drop_ka.nc")
                grid_l2(cfg_drop, spec, ts; footprint_method=SubpixelGridding(2),
                        outfile=out_drop_seq, backend=nothing, progress=false)
                grid_l2(cfg_drop, spec, ts; footprint_method=SubpixelGridding(2),
                        outfile=out_drop_ka, backend=ka_backend, progress=false)
                w_drop_seq = Dataset(out_drop_seq) do ds; Array(ds["n"][:, :, 1]); end
                w_drop_ka = Dataset(out_drop_ka) do ds; Array(ds["n"][:, :, 1]); end
                @test sum(w_drop_seq) == 0.0f0
                @test sum(w_drop_ka) == 0.0f0
                @test isequal(w_drop_ka, w_drop_seq)
            end

            # Case 4: Float32 fractional rows 1.0000004f0 / 5.0000005f0 with automatic oversampling (n=12)
            # Both direct fractional indices and end-to-end grid_l2, for quad and circle
            spec10 = GridSpec(lat_min=0f0, lat_max=10f0, lon_min=0f0, lon_max=10f0,
                              dlat=1f0, dlon=1f0)
            obs_auto = (; lat=Float32[2.0000004f0], lon=Float32[2.0f0],
                          lat_bounds=Float32[4.0f-7 4.0f-7 4.0000005f0 4.0000005f0],
                          lon_bounds=Float32[1.5f0 2.5f0 2.5f0 1.5f0])
            ilat_a, ilon_a, _ =
                SatelliteGridding._rectangular_indices(obs_auto, spec10, SubpixelGridding())
            @test ilat_a == Float32[1.0000004f0 1.0000004f0 5.0000005f0 5.0000005f0]
            auto_n_a = SatelliteGridding._l2_oversampling(ilat_a, ilon_a, nothing)
            @test auto_n_a == 12

            mktempdir() do dir
                f_auto = joinpath(dir, "auto_n_2020-01-01.nc")
                Dataset(f_auto, "c") do ds
                    defDim(ds, "x", 1)
                    defDim(ds, "nv", 4)
                    defVar(ds, "lat", Float32, ("x",))[:] = obs_auto.lat
                    defVar(ds, "lon", Float32, ("x",))[:] = obs_auto.lon
                    defVar(ds, "lat_bnds", Float32, ("x", "nv"))[:, :] = obs_auto.lat_bounds
                    defVar(ds, "lon_bnds", Float32, ("x", "nv"))[:, :] = obs_auto.lon_bounds
                    defVar(ds, "val", Float32, ("x",))[:] = Float32[13.0]
                end
                cfg_auto = DataSourceConfig(
                    Dict("lat" => "lat", "lon" => "lon",
                         "lat_bnd" => "lat_bnds", "lon_bnd" => "lon_bnds"),
                    OrderedDict("val" => "val"), FilterRule[],
                    "auto_n_YYYY-MM-DD.nc", dir)
                ts = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 1), Dates.Day(1))

                for (label, method) in (("quad", SubpixelGridding()),
                                        ("circle", CircularFootprintGridding()))
                    out_seq = joinpath(dir, "auto_$(label)_seq.nc")
                    out_ka = joinpath(dir, "auto_$(label)_ka.nc")
                    grid_l2(cfg_auto, spec10, ts; footprint_method=method,
                            outfile=out_seq, backend=nothing, progress=false)
                    grid_l2(cfg_auto, spec10, ts; footprint_method=method,
                            outfile=out_ka, backend=ka_backend, progress=false)
                    w_seq, v_seq = Dataset(out_seq) do ds
                        Array(ds["n"][:, :, 1]), Array(ds["val"][:, :, 1])
                    end
                    w_ka, v_ka = Dataset(out_ka) do ds
                        Array(ds["n"][:, :, 1]), Array(ds["val"][:, :, 1])
                    end
                    @test isequal(w_ka, w_seq)
                    mask = w_seq .> 0.0f0
                    @test any(mask)
                    @test all(isapprox.(v_ka[mask], v_seq[mask]; rtol=1e-5, atol=1e-5))
                end
            end

            # Case 5: L2 KA preserves degree-space CCW corner order instead of re-sorting
            # fractional indices (Float32 rounding can otherwise rotate the starting corner).
            mktempdir() do dir
                f_nosort = joinpath(dir, "nosort_2020-01-01.nc")
                Dataset(f_nosort, "c") do ds
                    defDim(ds, "x", 1)
                    defDim(ds, "nv", 4)
                    defVar(ds, "lat", Float32, ("x",))[:] = Float32[32.0f0]
                    defVar(ds, "lon", Float32, ("x",))[:] = Float32[10.0f0]
                    defVar(ds, "lat_bnds", Float32, ("x", "nv"))[:, :] =
                        Float32[prevfloat(32f0) 31.7f0 32.0f0 32.3f0]
                    defVar(ds, "lon_bnds", Float32, ("x", "nv"))[:, :] =
                        Float32[9.8f0 10.0f0 10.2f0 10.0f0]
                    defVar(ds, "val", Float32, ("x",))[:] = Float32[25.0f0]
                end
                cfg_nosort = DataSourceConfig(
                    Dict("lat" => "lat", "lon" => "lon",
                         "lat_bnd" => "lat_bnds", "lon_bnd" => "lon_bnds"),
                    OrderedDict("val" => "val"), FilterRule[],
                    "nosort_YYYY-MM-DD.nc", dir)
                ts = TimeSpec(DateTime(2020, 1, 1), DateTime(2020, 1, 1), Dates.Day(1))
                out_seq = joinpath(dir, "nosort_seq.nc")
                out_ka = joinpath(dir, "nosort_ka.nc")
                grid_l2(cfg_nosort, spec, ts; footprint_method=SubpixelGridding(2),
                        outfile=out_seq, backend=nothing, progress=false)
                grid_l2(cfg_nosort, spec, ts; footprint_method=SubpixelGridding(2),
                        outfile=out_ka, backend=ka_backend, progress=false)
                w_seq, v_seq = Dataset(out_seq) do ds
                    Array(ds["n"][:, :, 1]), Array(ds["val"][:, :, 1])
                end
                w_ka, v_ka = Dataset(out_ka) do ds
                    Array(ds["n"][:, :, 1]), Array(ds["val"][:, :, 1])
                end
                @test w_seq[96, 61:62] == Float32[0.25f0, 0.50f0]
                @test w_ka[96, 61:62] == Float32[0.25f0, 0.50f0]
                @test isequal(w_ka, w_seq)
                @test v_seq[96, 61:62] == Float32[25.0f0, 25.0f0]
                @test v_ka[96, 61:62] == Float32[25.0f0, 25.0f0]
            end
        end

        @testset "20,000 random soundings on 2° global grid — bounds, determinism, sequential match" begin
            rng = Test.Random.MersenneTwister(2026)
            spec = GridSpec(lat_min=-90f0, lat_max=90f0,
                            lon_min=-180f0, lon_max=180f0,
                            dlat=2f0, dlon=2f0)
            n_lon, n_lat = grid_shape(spec)  # (180, 90)
            N = 20_000
            lat = rand(rng, N) .* 180.0 .- 90.0
            lon = rand(rng, N) .* 360.0 .- 180.0
            # Scale lon offsets by 1/cos(lat) for realistic high-latitude footprint widths
            lon_scale = 1.0 ./ max.(abs.(cosd.(lat)), 1e-6)
            # Asymmetric corners (like real satellite footprints) filtered strictly inside domain
            lat_bnd = hcat(lat .- 0.3 .* rand(rng, N), lat .- 0.3 .* rand(rng, N),
                           lat .+ 0.3 .* rand(rng, N), lat .+ 0.3 .* rand(rng, N))
            lon_bnd = hcat(lon .- 0.3 .* rand(rng, N) .* lon_scale,
                           lon .+ 0.3 .* rand(rng, N) .* lon_scale,
                           lon .+ 0.3 .* rand(rng, N) .* lon_scale,
                           lon .- 0.3 .* rand(rng, N) .* lon_scale)
            sort_corners_ccw!(lat_bnd, lon_bnd)
            keep = (vec(minimum(lat_bnd, dims=2)) .> spec.lat_min) .&
                   (vec(maximum(lat_bnd, dims=2)) .< spec.lat_max) .&
                   (vec(minimum(lon_bnd, dims=2)) .> spec.lon_min) .&
                   (vec(maximum(lon_bnd, dims=2)) .< spec.lon_max) .&
                   ((vec(maximum(lon_bnd, dims=2)) .- vec(minimum(lon_bnd, dims=2))) .< 50.0)
            obs = (; lat=lat[keep], lon=lon[keep],
                   lat_bounds=lat_bnd[keep, :], lon_bounds=lon_bnd[keep, :])
            ilat64, ilon64, (ilat_c64, ilon_c64) =
                SatelliteGridding._rectangular_indices(obs, spec, CircularFootprintGridding())
            ilat32, ilon32, ilat_c32, ilon_c32 =
                map(a -> Float32.(a), (ilat64, ilon64, ilat_c64, ilon_c64))
            n_fp = length(ilat_c64)
            vals = randn(rng, Float32, n_fp, 2) .+ 10f0

            auto_n = SatelliteGridding._l2_oversampling(ilat64, ilon64, nothing)
            for n in (auto_n, 7)
                # 1. Circular index kernel produces no out-of-range indices
                ix = zeros(Int32, n_fp, n^2)
                iy = zeros(Int32, n_fp, n^2)
                ic = zeros(Int32, n_fp)
                sk = zeros(Int32, n_fp)
                compute_circular_footprint_indices_ka!(
                    ka_backend, ix, iy, ic, sk, ilat_c64, ilon_c64, ilat64, ilon64, n;
                    n_lon=n_lon, n_lat=n_lat)
                @test all(x -> x == Int32(-1) || (Int32(1) <= x <= Int32(n_lon)), ix)
                @test all(y -> y == Int32(-1) || (Int32(1) <= y <= Int32(n_lat)), iy)

                # 2. Run-to-run bit-determinism on KA CPU
                gs1 = zeros(Float32, n_lon, n_lat, 2)
                gw1 = zeros(Float32, n_lon, n_lat)
                accumulate_circular_batch!(ka_backend, gs1, gw1,
                                           ilat_c64, ilon_c64, ilat64, ilon64, vals, n, 2)

                gs2 = zeros(Float32, n_lon, n_lat, 2)
                gw2 = zeros(Float32, n_lon, n_lat)
                accumulate_circular_batch!(ka_backend, gs2, gw2,
                                           ilat_c64, ilon_c64, ilat64, ilon64, vals, n, 2)
                @test isequal(gw1, gw2)
                @test isequal(gs1, gs2)

                # 3. Per-cell circular agreement against sequential (Float64 & Float32 coords).
                # Populated-cell masks match identically; per-cell weights can differ where
                # KA's single-cell fast path adds 1.0f0 once vs summing rounded 1/inside_count.
                for (clat_in, clon_in, lat_in, lon_in) in (
                        (ilat_c64, ilon_c64, ilat64, ilon64),
                        (ilat_c32, ilon_c32, ilat32, ilon32))
                    gs_c = zeros(Float32, n_lon, n_lat, 2)
                    gw_c = zeros(Float32, n_lon, n_lat)
                    accumulate_circular_batch!(ka_backend, gs_c, gw_c,
                                               clat_in, clon_in, lat_in, lon_in, vals, n, 2)
                    gd_c_seq = zeros(Float32, n_lon, n_lat, 2)
                    gw_c_seq = zeros(Float32, n_lon, n_lat)
                    accumulate_circular_footprint!(gd_c_seq, nothing, gw_c_seq, false,
                                                   clat_in, clon_in, lat_in, lon_in, vals, n_fp, 2, n)
                    finalize_mean!(gs_c, gw_c)
                    @test (gw_c .> 0) == (gw_c_seq .> 0)
                    @test all(isapprox.(gw_c, gw_c_seq; rtol=1e-5, atol=2e-5))
                    mask = gw_c_seq .> 0.01f0
                    for z in 1:2
                        @test all(isapprox.(gs_c[:, :, z][mask], gd_c_seq[:, :, z][mask];
                                            rtol=1e-4, atol=1e-4))
                    end
                end

                # 4. Quadrilateral parity: exact per-cell weight equality (both Float64 & Float32 coords)
                pts = zeros(Float32, n, n, 2)
                ixb = zeros(Int32, n^2); iyb = zeros(Int32, n^2)
                l0 = zeros(n); o0 = zeros(n); l1 = zeros(n); o1 = zeros(n)
                for (lat_in, lon_in) in ((ilat64, ilon64), (ilat32, ilon32))
                    gs_q = zeros(Float32, n_lon, n_lat, 2)
                    gw_q = zeros(Float32, n_lon, n_lat)
                    accumulate_batch!(ka_backend, gs_q, gw_q,
                                      copy(lat_in), copy(lon_in), vals, n, 2)
                    gd_q_seq = zeros(Float32, n_lon, n_lat, 2)
                    gw_q_seq = zeros(Float32, n_lon, n_lat)
                    accumulate_footprint!(gd_q_seq, nothing, gw_q_seq, false,
                                          lat_in, lon_in, vals, n_fp, 2, n,
                                          pts, ixb, iyb, l0, o0, l1, o1)
                    finalize_mean!(gs_q, gw_q)
                    @test isequal(gw_q, gw_q_seq)
                    mask = gw_q_seq .> 0.01f0
                    for z in 1:2
                        @test all(isapprox.(gs_q[:, :, z][mask], gd_q_seq[:, :, z][mask];
                                            rtol=1e-4, atol=1e-4))
                    end
                end
            end
        end

        @testset "scatter_accumulate_circular! and scatter_accumulate_circular_ka! guard upper and lower bounds" begin
            n_lon, n_lat = 4, 4
            n = 2
            # fp 1 (flag=0): valid (2,2), col wrap OOB (5,1), row OOB (1,5), low OOB (0,1)
            # fp 2 (flag=1): col wrap OOB (5,1)
            # fp 3 (flag=1): row OOB (1,5)
            # fp 4 (flag=1): low OOB (0,2)
            # fp 5 (flag=1): valid (3,3)
            ix = Int32[2 5 1 0;
                       5 -1 -1 -1;
                       1 -1 -1 -1;
                       0 -1 -1 -1;
                       3 -1 -1 -1]
            iy = Int32[2 1 5 1;
                       1 -1 -1 -1;
                       5 -1 -1 -1;
                       2 -1 -1 -1;
                       3 -1 -1 -1]
            inside_count = Int32[4, 1, 1, 1, 1]
            skip = Int32[0, 1, 1, 1, 1]
            vals = Float32[12.0; 100.0; 200.0; 300.0; 7.0;;]

            # CPU sequential scatter
            gs_cpu = zeros(Float32, n_lon, n_lat, 1)
            gw_cpu = zeros(Float32, n_lon, n_lat)
            scatter_accumulate_circular!(gs_cpu, gw_cpu, ix, iy, inside_count, skip, vals, n, 1)

            # KA @atomic scatter (GPU path exercised on CPU backend)
            gs_ka = zeros(Float32, n_lon, n_lat, 1)
            gw_ka = zeros(Float32, n_lon, n_lat)
            scatter_accumulate_circular_ka!(ka_backend, gs_ka, gw_ka, ix, iy, inside_count, skip, vals, n, 1)

            for (gs, gw) in ((gs_cpu, gw_cpu), (gs_ka, gw_ka))
                @test gw[2, 2] ≈ 0.25f0
                @test gs[2, 2, 1] ≈ 3.0f0
                @test gw[3, 3] ≈ 1.0f0
                @test gs[3, 3, 1] ≈ 7.0f0
                # Specifically verify col=5,row=1 did NOT wrap into [1,2]
                @test gw[1, 2] == 0.0f0
                @test sum(gw) ≈ 1.25f0
                @test sum(gs) ≈ 10.0f0
            end
        end
    end

    @testset "Quadrilateral KA path — edge bounds and scatter guards" begin
        ka_backend = resolve_backend("cpu")

        @testset "compute_footprint_indices_ka! clips out-of-range indices" begin
            n_lon, n_lat = 180, 90
            # fp 1: fast path rounded to (181, 91) (OOB upper)
            # fp 2: fast path at (0, 1) (OOB lower)
            # fp 3: oversample straddling lat_max (90.5..91.5)
            # fp 4: oversample straddling lon_min (0.5..1.5)
            # fp 5: oversample fully in bounds
            lat_corners = Float32[
                91.0  91.0  91.0  91.0;
                1.2   1.2   1.8   1.8;
                90.5  90.5  91.5  91.5;
                45.0  45.0  46.0  46.0;
                10.0  10.0  12.0  12.0;
            ]
            lon_corners = Float32[
                181.0 181.0 181.0 181.0;
                0.2   0.8   0.8   0.2;
                90.0  91.0  91.0  90.0;
                0.5   1.5   1.5   0.5;
                20.0  22.0  22.0  20.0;
            ]
            n = 4
            ix = zeros(Int32, 5, n^2)
            iy = zeros(Int32, 5, n^2)
            skip = zeros(Int32, 5)

            compute_footprint_indices_ka!(ka_backend, ix, iy, skip, lat_corners, lon_corners, n;
                                          n_lon=n_lon, n_lat=n_lat)

            @test all(x -> x == Int32(-1) || (Int32(1) <= x <= Int32(n_lon)), ix)
            @test all(y -> y == Int32(-1) || (Int32(1) <= y <= Int32(n_lat)), iy)
            @test skip[1] == Int32(2)
            @test skip[2] == Int32(2)
            @test skip[3] == Int32(0)
            @test skip[4] == Int32(0)
            @test skip[5] == Int32(0)
        end

        @testset "scatter_accumulate! and scatter_accumulate_ka! guard upper and lower bounds" begin
            n_lon, n_lat = 4, 4
            n = 2
            ix = Int32[2 5 1 0;
                       5 -1 -1 -1;
                       1 -1 -1 -1;
                       0 -1 -1 -1;
                       3 -1 -1 -1]
            iy = Int32[2 1 5 1;
                       1 -1 -1 -1;
                       5 -1 -1 -1;
                       2 -1 -1 -1;
                       3 -1 -1 -1]
            skip = Int32[0, 1, 1, 1, 1]
            vals = Float32[12.0; 100.0; 200.0; 300.0; 7.0;;]

            gs_cpu = zeros(Float32, n_lon, n_lat, 1)
            gw_cpu = zeros(Float32, n_lon, n_lat)
            scatter_accumulate!(gs_cpu, gw_cpu, ix, iy, skip, vals, n, 1)

            gs_ka = zeros(Float32, n_lon, n_lat, 1)
            gw_ka = zeros(Float32, n_lon, n_lat)
            scatter_accumulate_ka!(ka_backend, gs_ka, gw_ka, ix, iy, skip, vals, n, 1)

            for (gs, gw) in ((gs_cpu, gw_cpu), (gs_ka, gw_ka))
                @test gw[2, 2] ≈ 0.25f0
                @test gs[2, 2, 1] ≈ 3.0f0
                @test gw[3, 3] ≈ 1.0f0
                @test gs[3, 3, 1] ≈ 7.0f0
                @test gw[1, 2] == 0.0f0
                @test sum(gw) ≈ 1.25f0
                @test sum(gs) ≈ 10.0f0
            end
        end

        @testset "accumulate_batch! drops out-of-grid subpixels cleanly" begin
            n_lon, n_lat = 6, 6
            n = 10
            # Straddles upper lat and lon edges, plus an out-of-bounds fast-path footprint
            lat_corners = Float32[
                5.5 5.5 7.5 7.5;
                7.0 7.0 7.0 7.0;
            ]
            lon_corners = Float32[
                5.5 7.5 7.5 5.5;
                7.0 7.0 7.0 7.0;
            ]
            vals = Float32[20.0; 99.0;;]

            gs1 = zeros(Float32, n_lon, n_lat, 1)
            gw1 = zeros(Float32, n_lon, n_lat)
            accumulate_batch!(ka_backend, gs1, gw1, copy(lat_corners), copy(lon_corners), vals, n, 1)

            gs2 = zeros(Float32, n_lon, n_lat, 1)
            gw2 = zeros(Float32, n_lon, n_lat)
            accumulate_batch!(ka_backend, gs2, gw2, copy(lat_corners), copy(lon_corners), vals, n, 1)

            @test isequal(gw1, gw2)
            @test isequal(gs1, gs2)

            expected_w = footprint_weights(SubpixelGridding(n),
                                           lat_corners[1, :], lon_corners[1, :];
                                           n_lon=n_lon, n_lat=n_lat)
            @test all(isapprox.(gw1, expected_w; rtol=1e-5, atol=1e-6))

            # Sequential accumulate_footprint! also guards bounds and matches KA bit-for-bit
            gd_seq = zeros(Float32, n_lon, n_lat, 1)
            gw_seq = zeros(Float32, n_lon, n_lat)
            pts = zeros(Float32, n, n, 2)
            ixb = zeros(Int32, n^2); iyb = zeros(Int32, n^2)
            l0 = zeros(n); o0 = zeros(n); l1 = zeros(n); o1 = zeros(n)
            accumulate_footprint!(gd_seq, nothing, gw_seq, false,
                                  lat_corners, lon_corners, vals, 2, 1, n,
                                  pts, ixb, iyb, l0, o0, l1, o1)
            @test isequal(gw1, gw_seq)

            finalize_mean!(gs1, gw1)
            mask = gw1 .> 0
            @test all(gs1[:, :, 1][mask] .≈ 20.0f0)
            @test all(gd_seq[:, :, 1][mask] .≈ 20.0f0)
        end
    end
end
