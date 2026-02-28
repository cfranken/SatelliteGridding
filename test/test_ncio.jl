@testset "NetCDF I/O" begin

    @testset "read_nc_variable — flat variable" begin
        tmpfile = tempname() * ".nc"
        ds = Dataset(tmpfile, "c")
        defDim(ds, "x", 5)
        v = defVar(ds, "temperature", Float32, ("x",))
        v[:] = Float32[1.0, 2.0, 3.0, 4.0, 5.0]
        close(ds)

        ds = Dataset(tmpfile)
        data = read_nc_variable(ds, "temperature")
        @test data ≈ Float32[1.0, 2.0, 3.0, 4.0, 5.0]
        close(ds)
        rm(tmpfile)
    end

    @testset "read_nc_variable — nested groups" begin
        tmpfile = tempname() * ".nc"
        ds = Dataset(tmpfile, "c")
        defDim(ds, "x", 3)
        grp = defGroup(ds, "PRODUCT")
        v = defVar(grp, "methane", Float32, ("x",))
        v[:] = Float32[1800.0, 1850.0, 1900.0]
        close(ds)

        ds = Dataset(tmpfile)
        data = read_nc_variable(ds, "PRODUCT/methane")
        @test data ≈ Float32[1800.0, 1850.0, 1900.0]
        close(ds)
        rm(tmpfile)
    end

    @testset "read_nc_variable — bounds reshaping (4 first)" begin
        tmpfile = tempname() * ".nc"
        ds = Dataset(tmpfile, "c")
        defDim(ds, "nv", 4)
        defDim(ds, "x", 3)
        v = defVar(ds, "lat_bnds", Float32, ("nv", "x"))
        v[:] = Float32[1 4 7; 2 5 8; 3 6 9; 4 7 10]  # 4×3
        close(ds)

        ds = Dataset(tmpfile)
        data = read_nc_variable(ds, "lat_bnds"; bounds=true)
        @test size(data) == (3, 4)  # N×4
        @test data[1, :] ≈ Float32[1, 2, 3, 4]
        close(ds)
        rm(tmpfile)
    end

    @testset "read_nc_variable — bounds reshaping (4 last)" begin
        tmpfile = tempname() * ".nc"
        ds = Dataset(tmpfile, "c")
        defDim(ds, "x", 3)
        defDim(ds, "nv", 4)
        v = defVar(ds, "lon_bnds", Float32, ("x", "nv"))
        v[:] = Float32[1 2 3 4; 5 6 7 8; 9 10 11 12]  # 3×4
        close(ds)

        ds = Dataset(tmpfile)
        data = read_nc_variable(ds, "lon_bnds"; bounds=true)
        @test size(data) == (3, 4)
        @test data[1, :] ≈ Float32[1, 2, 3, 4]
        close(ds)
        rm(tmpfile)
    end

    @testset "read_nc_attribute" begin
        tmpfile = tempname() * ".nc"
        ds = Dataset(tmpfile, "c")
        defDim(ds, "x", 3)
        v = defVar(ds, "sif", Float32, ("x",),
                   attrib=["units" => "mW/m²/sr/nm", "long_name" => "Solar Induced Fluorescence"])
        v[:] = Float32[0.5, 1.0, 1.5]
        close(ds)

        ds = Dataset(tmpfile)
        @test read_nc_attribute(ds, "sif", "units") == "mW/m²/sr/nm"
        @test read_nc_attribute(ds, "sif", "long_name") == "Solar Induced Fluorescence"
        close(ds)
        rm(tmpfile)
    end

    @testset "create_output_dataset" begin
        tmpfile = tempname() * ".nc"
        gs = GridSpec(lat_min=-5.0f0, lat_max=5.0f0,
                      lon_min=-10.0f0, lon_max=10.0f0,
                      dlat=5.0f0, dlon=5.0f0)
        grid_vars = OrderedDict("sif" => "SIF_757nm", "cloud" => "cloud_fraction")

        ds, nc_vars = create_output_dataset(tmpfile, gs, 3, grid_vars, true)

        @test haskey(nc_vars, "sif")
        @test haskey(nc_vars, "sif_std")
        @test haskey(nc_vars, "cloud")
        @test haskey(nc_vars, "cloud_std")
        @test haskey(nc_vars, "n")
        @test size(ds["lat"]) == (2,)
        @test size(ds["lon"]) == (4,)
        @test size(ds["time"]) == (3,)

        close(ds)
        rm(tmpfile)
    end
end
