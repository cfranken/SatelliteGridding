@testset "NetCDF I/O" begin

    @testset "Bounds have one canonical orientation, including four soundings" begin
        mktempdir() do dir
            for n in (3, 4, 5), vertex_first in (false, true), grouped in (false, true)
                file = joinpath(dir, "bounds.nc")
                expected = Float32[10i + corner for i in 1:n, corner in (0.1, 0.1, 0.9, 0.9)]
                Dataset(file, "c") do ds
                    defDim(ds, "sounding", n)
                    defDim(ds, "vertex", 4)
                    parent = grouped ? defGroup(ds, "Geolocation") : ds
                    dims = vertex_first ? ("vertex", "sounding") : ("sounding", "vertex")
                    bounds = defVar(parent, "latitude_bounds", Float32, dims)
                    bounds[:, :] = vertex_first ? permutedims(expected) : expected
                end
                Dataset(file) do ds
                    path = grouped ? "Geolocation/latitude_bounds" : "latitude_bounds"
                    @test read_nc_variable(ds, path; bounds=true) == expected
                    @test size(read_nc_variable(ds, path; bounds=true)) == (n, 4)
                end
            end

            # Singleton time dimensions and an interior vertex axis also flatten
            # into sounding × vertex, using the same path as ordinary bounds.
            file = joinpath(dir, "nested.nc")
            expected = Float32[10i + corner for i in 1:4, corner in 1:4]
            Dataset(file, "c") do ds
                defDim(ds, "time", 1)
                defDim(ds, "sounding", 4)
                defDim(ds, "vertex", 4)
                group = defGroup(ds, "Geo")
                v = defVar(group, "bounds", Float32, ("time", "vertex", "sounding"))
                v[:, :, :] = reshape(permutedims(expected), 1, 4, 4)
            end
            Dataset(file) do ds
                @test read_nc_variable(ds, "Geo/bounds"; bounds=true) == expected
            end
        end
    end

    @testset "Output definition failures close the dataset immediately" begin
        mktempdir() do dir
            opened = Ref{Any}(nothing)
            @test_throws ErrorException SatelliteGridding._with_output_setup(joinpath(dir, "failed.nc")) do ds
                opened[] = ds
                defDim(ds, "x", 1)
                error("synthetic failure while defining output")
            end
            @test !NCDatasets.isopen(opened[])
            # On success, the caller owns the open dataset, as before.
            ds = SatelliteGridding._with_output_setup(identity, joinpath(dir, "ok.nc"))
            @test NCDatasets.isopen(ds)
            close(ds)
        end
    end

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

        filters = [FilterRule("qa_value", :gt, 0.5),
                   FilterRule("methane", :between, 1600.0, 2200.0)]
        ds, nc_vars = create_output_dataset(tmpfile, gs, 3, grid_vars, true;
                                            command="julia bin/grid.jl l2 --config data.toml",
                                            filters=filters)

        @test haskey(nc_vars, "sif")
        @test haskey(nc_vars, "sif_std")
        @test haskey(nc_vars, "cloud")
        @test haskey(nc_vars, "cloud_std")
        @test haskey(nc_vars, "n")
        @test size(ds["lat"]) == (2,)
        @test size(ds["lon"]) == (4,)
        @test size(ds["time"]) == (3,)
        @test ds.attrib["command"] == "julia bin/grid.jl l2 --config data.toml"
        @test ds.attrib["software_version"] == "0.1.0"
        @test ds.attrib["julia_version"] == string(VERSION)
        @test haskey(ds.attrib, "git_tag")
        @test haskey(ds.attrib, "git_commit")
        @test haskey(ds.attrib, "git_dirty")
        @test ds.attrib["filters"] ==
              "qa_value: x > 0.5; methane: 1600.0 < x < 2200.0"

        close(ds)
        rm(tmpfile)
    end

    @testset "NetCDF storage settings survive reopening" begin
        mktempdir() do dir
            gs = GridSpec(lat_min=0f0, lat_max=2f0, lon_min=0f0, lon_max=3f0, dlat=1f0, dlon=1f0)
            gv = OrderedDict("value" => "value")
            for (name, options, level, shuffled) in (
                ("default", NetCDFOptions(), 4, false),
                ("off", NetCDFOptions(compression="none"), 0, false),
                ("fast", NetCDFOptions(deflate_level=1, shuffle=false, fill_value=-1234), 1, false),
                ("small", NetCDFOptions(deflate_level=9, shuffle=true, fill_value=NaN, chunksizes=[2, 2, 1],
                                         checksum="fletcher32"), 9, true),
            )
                file = joinpath(dir, name * ".nc")
                ds, vars = create_output_dataset(file, gs, 2, gv, true; output_options=options)
                # Include populated, missing, and entirely unwritten time slices.
                vars["value"][:, :, 1] = Float32[1 2; options.fill_value 4; 5 6]
                close(ds)
                Dataset(file) do reopened
                    for key in ("value", "value_std", "n")
                        v = reopened[key]
                        shuf, deflated, lev = NCDatasets.deflate(v.var)
                        @test lev == level
                        @test shuf == shuffled
                        @test deflated == (level > 0)
                        @test isequal(v.attrib["_FillValue"], options.fill_value)
                        expected_chunks = options.chunksizes === nothing ? (3, 2, 1) : options.chunksizes
                        @test NCDatasets.chunking(v.var)[2] == expected_chunks
                        @test NCDatasets.checksum(v.var) == options.checksum
                    end
                    @test reopened["value"][1, 1, 1] == 1
                    @test ismissing(reopened["value"][2, 1, 1])
                    @test all(ismissing, reopened["value"][:, :, 2])
                end
            end

            cube = CubedSphereGridSpec(Nc=4)
            file = joinpath(dir, "cube.nc")
            ds, vars = create_output_dataset(file, cube, 1, gv, true;
                output_options=NetCDFOptions(fill_value=NaN, chunksizes=[2, 2, 1, 1]))
            @test NCDatasets.chunking(vars["value"].var)[2] == (2, 2, 1, 1)
            @test isnan(vars["value_std"].attrib["_FillValue"])
            close(ds)
            @test_throws ArgumentError create_output_dataset(joinpath(dir, "invalid.nc"), gs, 1, gv, false;
                output_options=NetCDFOptions(chunksizes=[3, 2]))
            @test !isfile(joinpath(dir, "invalid.nc"))
        end
    end

    @testset "Git provenance, including source archives" begin
        mktempdir() do dir
            run(pipeline(`git -C $dir init`; stdout=devnull, stderr=devnull))
            write(joinpath(dir, "code.txt"), "original\n")
            run(`git -C $dir add code.txt`)
            run(pipeline(`git -C $dir -c user.name=Test -c user.email=test@example.invalid commit -m initial`;
                         stdout=devnull, stderr=devnull))
            run(`git -C $dir tag v1.2.3`)
            provenance = SatelliteGridding._software_provenance(; root=dir)
            @test provenance["git_tag"] == "v1.2.3"
            @test length(provenance["git_commit"]) == 40
            @test provenance["git_dirty"] == "false"
            write(joinpath(dir, "code.txt"), "changed\n")
            @test SatelliteGridding._software_provenance(; root=dir)["git_dirty"] == "true"
        end
        mktempdir() do dir
            provenance = SatelliteGridding._software_provenance(; root=dir)
            @test provenance["git_tag"] == "unknown"
            @test provenance["git_commit"] == "unknown"
            @test provenance["git_dirty"] == "unknown"
        end
    end

    @testset "OCO version layouts and metadata selection" begin
        mktempdir() do dir
            file = joinpath(dir, "metadata.nc")
            Dataset(file, "c") do ds
                ds.attrib["BuildId"] = "B11.3.00"
                ds.attrib["institution"] = "JPL"
                group = defGroup(ds, "Metadata")
                group.attrib["version"] = "group-version"
                defVar(group, "BuildId", String, ())[] = "B11.2.17"
            end
            cfg = DataSourceConfig(Dict{String,String}(), OrderedDict("a" => "a"),
                                    FilterRule[], "*.nc", dir)
            metadata = SatelliteGridding.InputMetadata(cfg)
            Dataset(file) do ds
                SatelliteGridding._capture_input_metadata!(metadata, ds)
            end
            @test metadata.product_versions == Set(["B11.3.00"])
            @test metadata.versions["Metadata/BuildId"] == ["B11.2.17"]
            cfg.options["metadata"] = Dict("version_sources" => ["Metadata/@version"],
                "attributes" => ["@institution", "Metadata/BuildId", "@missing"])
            selected = SatelliteGridding.InputMetadata(cfg)
            Dataset(file) do ds
                SatelliteGridding._capture_input_metadata!(selected, ds)
                SatelliteGridding._capture_input_metadata!(selected, ds)
            end
            @test selected.product_versions == Set(["group-version"])
            @test selected.attributes["@institution"] == ["JPL"]
            @test selected.attributes["Metadata/BuildId"] == ["B11.2.17"]
            @test !haskey(selected.attributes, "@missing")
        end
    end
end
