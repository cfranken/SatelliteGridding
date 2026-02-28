@testset "Filters" begin

    @testset "apply_center_filters" begin
        lat = Float32[5 15 25; 35 45 55]
        lon = Float32[10 20 30; 40 50 60]

        gs = GridSpec(lat_min=10.0f0, lat_max=50.0f0,
                      lon_min=10.0f0, lon_max=50.0f0,
                      dlat=10.0f0, dlon=10.0f0)

        idx = SatelliteGridding.apply_center_filters(lat, lon, gs)

        # Only pixels strictly within (10,50) x (10,50):
        # (15,20), (25,30), (35,40), (45,50) — but 50 is not < 50, so (45,50) excluded
        # Check that boundary pixels are excluded
        for ci in idx
            @test lat[ci] > 10.0f0
            @test lat[ci] < 50.0f0
            @test lon[ci] > 10.0f0
            @test lon[ci] < 50.0f0
        end
    end
end
