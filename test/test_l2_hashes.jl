using SHA

function _write_l2_hash_fixture(dir)
    rectangles = [(32.1, 49.3, -2.2, 2.6), (33.2, 50.0, -1.4, 2.2),
                  (0.1, 1.1, 0.75, 1.25), (0.5, 1.6, 0.25, 1.75),
                  (359.5, 360.8, 0.15, 1.45), (321.1, 322.3, 88.8, 90.0)]
    for day in (1, 2, 3, 5)
        Dataset(joinpath(dir, "input_2020-01-0$day.nc"), "c") do ds
            ds.attrib["product_version"] = "fixture-v$day"
            defDim(ds, "obs", length(rectangles))
            defDim(ds, "vertex", 4)
            lat = defVar(ds, "lat", Float64, ("obs", "vertex"))
            lon = defVar(ds, "lon", Float64, ("obs", "vertex"))
            clat = defVar(ds, "clat", Float64, ("obs",))
            clon = defVar(ds, "clon", Float64, ("obs",))
            for (i, (x0, x1, y0, y1)) in enumerate(rectangles)
                dx = 0.11 * day
                lat[i, :] = [y1, y0, y1, y0]
                lon[i, :] = [x0, x1, x1, x0] .+ dx
                clat[i] = (y0 + y1) / 2
                clon[i] = (x0 + x1) / 2 + dx
            end
            defVar(ds, "radius", Float64, ("obs",))[:] = [0.8, 1.1, 0.5, 0.71, 1.2, 1.0]
            defVar(ds, "signal", Float64, ("obs",); attrib=["units" => "K"])[:] =
                [1.123456789, -4.3456789, 13.56789123, 31.78912345, 3.141592653, 0.023456789] .+ (day / 7)
            defVar(ds, "offset", Float64, ("obs",))[:] = 10_000.0 .+ ((1:6) .+ day) ./ 989.0
            defVar(ds, "qa", Float64, ("obs",))[:] =
                day == 5 ? zeros(6) : [0.99, 0.9, 1.0, 0.89, 0.98, 0.79]
        end
    end
end

function _write_l2_hash_data!(io, data::AbstractArray{T}) where {T}
    print(io, T, size(data))
    if T <: Union{Float32,Float64,Int32,Int64,UInt32,UInt64}
        U = sizeof(T) == 4 ? UInt32 : UInt64
        write(io, reinterpret(UInt8, htol.(reinterpret(U, vec(data)))))
    elseif T === Char
        write(io, string(vec(data)...))
    else
        for value in data
            print(io, repr(value), '\n')
        end
    end
end
_write_l2_hash_data!(io, ::Nothing) = print(io, "nothing")
_write_l2_hash_data!(io, value) = print(io, repr(value))

function _l2_payload_hash(path)
    payload = IOBuffer()
    Dataset(path) do ds
        for key in keys(ds)
            print(payload, key)
            raw = ds[key].var
            _write_l2_hash_data!(payload, ndims(raw) == 0 ? raw[] : raw[ntuple(_ -> Colon(), ndims(raw))...])
        end
    end
    bytes2hex(sha256(take!(payload)))
end

@testset "L2 NetCDF payload bit regressions" begin
    # Raw variable hashes from pre-cleanup HEAD, 11115f6. The fixture
    # includes a QA-rejected file, empty days, overlapping footprints, and
    # polar and seam geometry. Hash the on-disk values without CF decoding.
    cases = (
        ("rect", Float32, "quad", false, true, "float32") => "6a1ff8b417621839f23276c05541ff190839722c6feab688ec7c2b8fbf88cf04",
        ("rect", Float64, "circle", true, true, "float64") => "62ed7d54bc77c6eef174c934f7bd85d29061512300299f604a879936e7eabcec",
        ("rect", Float64, "radius", false, false, "float32") => "2c95678b9e52c862e0d6c2752d45db5efa0f81ab40ace4e6f94d22a9df320322",
        ("rect-ka", Float64, "quad", false, false, "float32") => "713a9680afd94f678d0baee28e4e939367fdd77104940184631c7f12d99f1897",
        ("rect-ka", Float32, "quad", true, false, "float64") => "3b1b5a007d40ba9195d72128dd1bef2e1ea24623bac662a1c8153b229b8b7ef9",
        ("rect-ka", Float32, "circle", true, false, "float64") => "d072c13c2a1b29ef3cbbc5adfd561d368587c7fc1f7aa5871c255173a83fc234",
        ("rect-ka", Float64, "radius", false, false, "float64") => "dcb11a185ce9ada5d9fdaee696baa666cc4b3871ff877ea8586d5e09922275e2",
        ("cube", Float32, "quad", true, true, "float32") => "f154a166e135e85e0a726d538d9f5b9d802498a2d16346e4d83b509541857c76",
        ("cube-eq", Float64, "quad", false, false, "float64") => "f6cfca012ef16a7f27703b268fc27a367c8256c480f41afe13322b0e8cf092dc",
        ("healpix", Float32, "quad", false, true, "float64") => "060b7db41885239bd28fa90ef6237011766657c890d740e8cba8298427ab46bf",
        ("healpix", Float64, "quad", true, false, "float32") => "aa03ff9e3216e0baebe87b0c75bbbef4461e9772c985ef2fe040bb663700174b",
        ("healpix-ka", Float32, "quad", true, true, "float32") => "00033007e3941cbeeb97e9b60a9cc60f5b20bcd4409fd672457585ad3e15eede",
        ("healpix-ka", Float64, "quad", false, true, "float64") => "ffe6e0cbe638e9d0ff30cfefdf73a5d0cb542cc689a353e1846fc16ecf3bbb17",
    )
    scratch = mkpath(joinpath(@__DIR__, "..", ".scratch"))
    mktempdir(scratch) do dir
        _write_l2_hash_fixture(dir)
        basic = Dict("lat_bnd" => "lat", "lon_bnd" => "lon", "lat" => "clat",
                     "lon" => "clon", "radius" => "radius")
        vars = OrderedDict("signal" => "signal", "offset" => "offset")
        filters = [FilterRule("qa", :gt, 0.8)]
        cfg = DataSourceConfig(basic, vars, filters, "input_YYYY-MM-DD.nc", dir)
        radius_cfg = DataSourceConfig(Dict(k => v for (k, v) in basic if !endswith(k, "_bnd")),
                                      vars, filters, cfg.file_pattern, dir)
        start, stop = DateTime(2020, 1, 1), DateTime(2020, 1, 6)
        tiled = TimeSpec(start, stop, Day(1); oversample_temporal=2.0f0)
        rolling = TimeSpec(start, stop, Day(1); sample_step=Day(1), window_halfwidth_days=1)
        path = joinpath(dir, "output.nc")
        for ((label, T, footprint, roll, std, output), expected) in cases
            @testset "$label $T $footprint rolling=$roll std=$std output=$output" begin
                spec = if startswith(label, "rect")
                    GridSpec(lat_min=T(-5), lat_max=T(5), lon_min=T(0), lon_max=T(60), dlat=T(1), dlon=T(1))
                elseif startswith(label, "cube")
                    def = label == "cube" ? GMAOCubedSphereDefinition() : EquiangularCubedSphereDefinition()
                    CubedSphereGridSpec(Nc=12, T=T, definition=def)
                else
                    HEALPixGridSpec(nside=8, T=T)
                end
                config = footprint == "radius" ? radius_cfg : cfg
                config.options["output"] = output == "float32" ? Dict("precision" => output) :
                    Dict("precision" => output, "significant_digits_mean" => "none", "significant_digits_std" => "none")
                method = footprint == "quad" ? SubpixelGridding() : CircularFootprintGridding(7)
                backend = endswith(label, "-ka") ? CPU() : nothing
                grid_l2(config, spec, roll ? rolling : tiled; n_oversample=7,
                        footprint_method=method, outfile=path, backend, compute_std=std,
                        progress=false, verbose=false)
                @test _l2_payload_hash(path) == expected
            end
        end

        @testset "Rectangular KA STD fails before opening output" begin
            sentinel = read(path)
            spec = GridSpec(lat_min=-5f0, lat_max=5f0, lon_min=0f0, lon_max=60f0, dlat=1f0, dlon=1f0)
            message = "Rectangular KA gridding does not support compute_std=true; use --backend sequential or omit --compSTD."
            for schedule in (tiled, rolling)
                @test_throws ErrorException(message) grid_l2(
                    cfg, spec, schedule; outfile=path, backend=CPU(), compute_std=true, progress=false, verbose=false)
                @test read(path) == sentinel
            end
        end

        @testset "An empty selection still counts as a successful file" begin
            write(joinpath(dir, "input_2020-01-04.nc"), "deliberately invalid NetCDF")
            day4, day5 = DateTime(2020, 1, 4), DateTime(2020, 1, 5)
            schedules = (TimeSpec(day4, day4, Day(2)),
                         TimeSpec(day5, day5, Day(1); sample_step=Day(1), window_halfwidth_days=1))
            for schedule in schedules
                @test_logs (:warn, "Error processing file; continuing") grid_l2(
                    cfg, HEALPixGridSpec(nside=8), schedule; outfile=path, backend=CPU(),
                    compute_std=true, n_oversample=7, keep_going=true, progress=false, verbose=false)
                Dataset(path) do ds
                    @test all(iszero, ds["n"][:])
                    @test all(ismissing, ds["signal_std"][:])
                end
            end
        end
    end
end
