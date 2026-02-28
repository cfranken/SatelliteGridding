@testset "Oversampling" begin

    @testset "subdivide_line!" begin
        # Line from (0,0) to (10,10) with n=5
        # Sub-segments have width 2, centers at 1, 3, 5, 7, 9
        points = zeros(Float32, 1, 5, 2)
        SatelliteGridding.subdivide_line!(points, 1, 0.0, 0.0, 10.0, 10.0, 5)
        @test points[1, :, 1] ≈ [1.0, 3.0, 5.0, 7.0, 9.0]
        @test points[1, :, 2] ≈ [1.0, 3.0, 5.0, 7.0, 9.0]

        # Line from (0,0) to (0,6) with n=3 (only lon changes)
        points2 = zeros(Float32, 1, 3, 2)
        SatelliteGridding.subdivide_line!(points2, 1, 0.0, 0.0, 0.0, 6.0, 3)
        @test all(points2[1, :, 1] .≈ 0.0)  # lat stays 0
        @test points2[1, :, 2] ≈ [1.0, 3.0, 5.0]
    end

    @testset "subdivide_baseline!" begin
        lats = zeros(4)
        lons = zeros(4)
        SatelliteGridding.subdivide_baseline!(lats, lons, 0.0, 0.0, 8.0, 8.0, 4)
        @test lats ≈ [1.0, 3.0, 5.0, 7.0]
        @test lons ≈ [1.0, 3.0, 5.0, 7.0]
    end

    @testset "compute_subpixels! — unit square" begin
        # Footprint: unit square (0,0), (0,1), (1,1), (1,0)
        # Corner order: 1=(0,0), 2=(0,1), 3=(1,1), 4=(1,0)
        n = 3
        points = zeros(Float32, n, n, 2)
        lats0 = zeros(n)
        lons0 = zeros(n)
        lats1 = zeros(n)
        lons1 = zeros(n)

        vert_lat = Float32[0, 0, 1, 1]
        vert_lon = Float32[0, 1, 1, 0]

        compute_subpixels!(points, vert_lat, vert_lon, n, lats0, lons0, lats1, lons1)

        # With n=3 on a unit square, sub-pixels should be at 1/6, 3/6, 5/6 in each dimension
        # First index (i) iterates along edge 1→2 (lon direction), second (j) along connecting lines (lat direction)
        expected = [1/6, 3/6, 5/6]
        # lat varies along j (second dim), constant along i (first dim)
        for j in 1:3
            @test points[:, j, 1] ≈ fill(Float32(expected[j]), 3) atol=1e-5
        end
        # lon varies along i (first dim), constant along j (second dim)
        for i in 1:3
            @test points[i, :, 2] ≈ fill(Float32(expected[i]), 3) atol=1e-5
        end

        # All points should be inside the unit square
        @test all(0 .< points[:, :, 1] .< 1)
        @test all(0 .< points[:, :, 2] .< 1)
    end

    @testset "compute_subpixels! — rectangular footprint" begin
        # 2° wide, 1° tall footprint
        n = 4
        points = zeros(Float32, n, n, 2)
        lats0 = zeros(n); lons0 = zeros(n)
        lats1 = zeros(n); lons1 = zeros(n)

        vert_lat = Float32[0, 0, 1, 1]
        vert_lon = Float32[0, 2, 2, 0]

        compute_subpixels!(points, vert_lat, vert_lon, n, lats0, lons0, lats1, lons1)

        # All 16 sub-pixels within footprint bounds
        @test all(0 .< points[:, :, 1] .< 1)
        @test all(0 .< points[:, :, 2] .< 2)
    end

    @testset "compute_subpixels! — tilted quadrilateral" begin
        # Non-axis-aligned diamond shape
        n = 5
        points = zeros(Float32, n, n, 2)
        lats0 = zeros(n); lons0 = zeros(n)
        lats1 = zeros(n); lons1 = zeros(n)

        vert_lat = Float32[0, 1, 2, 1]
        vert_lon = Float32[1, 0, 1, 2]

        compute_subpixels!(points, vert_lat, vert_lon, n, lats0, lons0, lats1, lons1)

        # Sub-pixels should be within the bounding box of the diamond
        @test all(-0.1 .< points[:, :, 1] .< 2.1)
        @test all(-0.1 .< points[:, :, 2] .< 2.1)
    end

    @testset "floor_indices!" begin
        points = Float32[1.7 2.3; 3.1 4.9]
        points3d = reshape(points, 1, 2, 2)  # (1, 2, 2)
        ix = zeros(Int32, 2)
        iy = zeros(Int32, 2)
        floor_indices!(ix, iy, points3d)
        @test ix == Int32[1, 3]
        @test iy == Int32[2, 4]
    end

    @testset "compute_n_oversample" begin
        # Small footprint relative to grid → small n
        @test compute_n_oversample(0.5, 1.0) == 2  # clamped to n_min

        # Footprint ≈ grid cell → moderate n
        @test compute_n_oversample(1.0, 1.0) == 3

        # Large footprint → large n
        @test compute_n_oversample(5.0, 1.0) == 15

        # Very large footprint → clamped to n_max
        @test compute_n_oversample(20.0, 1.0) == 20

        # Custom bounds
        @test compute_n_oversample(0.1, 1.0; n_min=5) == 5
    end

    @testset "sort_corners_ccw" begin
        # Already counter-clockwise (like TROPOMI)
        lat_ccw = Float32[-1, -1, 1, 1]
        lon_ccw = Float32[-1, 1, 1, -1]
        lat_s, lon_s = sort_corners_ccw(lat_ccw, lon_ccw)
        # Should maintain CCW order (angles: SW≈-135, SE≈-45, NE≈45, NW≈135)
        @test lat_s ≈ Float32[-1, -1, 1, 1]
        @test lon_s ≈ Float32[-1, 1, 1, -1]

        # Scrambled order (like OCO-2: 1,4,3,2 → SW,NW,NE,SE)
        lat_scr = Float32[-1, 1, 1, -1]
        lon_scr = Float32[-1, -1, 1, 1]
        lat_s2, lon_s2 = sort_corners_ccw(lat_scr, lon_scr)
        # After sorting, should be CCW: SW(-1,-1), SE(-1,1), NE(1,1), NW(1,-1)
        @test lat_s2 ≈ Float32[-1, -1, 1, 1]
        @test lon_s2 ≈ Float32[-1, 1, 1, -1]
    end

    @testset "sort_corners_ccw!" begin
        # Test in-place version with N×4 matrix
        lat_mat = Float32[-1 1 1 -1; 0 1 0 -1]
        lon_mat = Float32[-1 -1 1 1; -1 0 1 0]
        sort_corners_ccw!(lat_mat, lon_mat)

        # First sounding: scrambled → CCW
        @test lat_mat[1, :] ≈ Float32[-1, -1, 1, 1]
        @test lon_mat[1, :] ≈ Float32[-1, 1, 1, -1]

        # Second sounding: diamond → CCW (bottom, right, top, left)
        @test lat_mat[2, :] ≈ Float32[-1, 0, 1, 0]
        @test lon_mat[2, :] ≈ Float32[0, 1, 0, -1]
    end
end
