@testset "Filters" begin

    @testset "Packed quality flags and missing values" begin
        mktempdir() do dir
            file = joinpath(dir, "quality.nc")
            Dataset(file, "c") do ds
                defDim(ds, "x", 4)
                qa = defVar(ds, "qa", Int16, ("x",); fillvalue=Int16(32767),
                            attrib=["scale_factor" => 0.01f0])
                qa.var[:] = Int16[40, 60, 32767, 90]
            end
            cfg = DataSourceConfig(Dict{String,String}(), OrderedDict("a" => "a"),
                [FilterRule("qa", :gt, 0.5)], "*.nc", dir)
            lat_bounds = repeat(Float32[0.25 0.25 0.75 0.75], 4, 1)
            lon_bounds = repeat(Float32[0.25 0.75 0.75 0.25], 4, 1)
            gs = GridSpec(lat_min=0f0, lat_max=1f0, lon_min=0f0, lon_max=1f0, dlat=1f0, dlon=1f0)
            Dataset(file) do ds
                @test apply_filters(ds, cfg, lat_bounds, lon_bounds, gs) == [2, 4]
                @test apply_filters(ds, cfg, lat_bounds, lon_bounds, CubedSphereGridSpec(Nc=2)) == [2, 4]
            end
        end
    end

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
