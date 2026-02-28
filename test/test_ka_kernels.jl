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

        @test lat_ka ≈ lat_mat
        @test lon_ka ≈ lon_mat
    end

    @testset "sort_corners_ccw_ka! — many soundings" begin
        n = 1000
        lat_mat = randn(Float32, n, 4)
        lon_mat = randn(Float32, n, 4)

        lat_ka = copy(lat_mat)
        lon_ka = copy(lon_mat)

        sort_corners_ccw!(lat_mat, lon_mat)
        sort_corners_ccw_ka!(backend, lat_ka, lon_ka)

        @test lat_ka ≈ lat_mat
        @test lon_ka ≈ lon_mat
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

        # Compare weights
        @test grid_w_ka ≈ grid_w_seq rtol=1e-4

        # Compare means where weight > 0
        mask = grid_w_seq .> 0.01f0
        for z in 1:n_vars
            ka_vals = grid_sum[:, :, z][mask]
            seq_vals = grid_data[:, :, z][mask]
            @test ka_vals ≈ seq_vals rtol=1e-3
        end
    end

    @testset "finalize_mean!" begin
        grid_sum = Float32[10.0;;; 20.0]  # shape (1,1,2)
        grid_weights = Float32[5.0;;]     # shape (1,1)

        finalize_mean!(grid_sum, grid_weights)

        @test grid_sum[1, 1, 1] ≈ 2.0f0   # 10/5
        @test grid_sum[1, 1, 2] ≈ 4.0f0   # 20/5
    end
end
