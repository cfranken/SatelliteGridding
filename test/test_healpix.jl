using Random
using SparseArrays

_healpix_ref_rows(filename) =
    (split(line, ',') for line in Iterators.drop(readlines(joinpath(@__DIR__, "data", filename)), 1))

_healpix_unit_vector(lon, lat) = (cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))
function _healpix_angle_error(lon, lat, lon_ref, lat_ref)
    a = _healpix_unit_vector(lon, lat)
    b = _healpix_unit_vector(lon_ref, lat_ref)
    chord = sqrt(sum(abs2, a .- b))
    return 2 * asind(chord / 2)
end

# Test footprints are supplied in perimeter order. Evaluate the four direct
# bilinear basis functions independently of production's corner sort and its
# two-stage interpolation. Angles derived from component ratios are scale-free.
function _healpix_reference_counts(nside, lon, lat, n; spherical=true)
    vertices = spherical ? map(_healpix_unit_vector, lon, lat) : map(tuple, lon, lat)
    counts = Dict{Int,Int}()
    for j in 1:n, i in 1:n
        u, v = (i - 0.5) / n, (j - 0.5) / n
        basis = ((1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v)
        x, y = ntuple(axis -> sum(basis[k] * vertices[k][axis] for k in 1:4), 2)
        if spherical
            z = sum(basis[k] * vertices[k][3] for k in 1:4)
            x, y = rad2deg(atan(y, x)), rad2deg(atan(z, hypot(x, y)))
        end
        pix = ang2pix_nest(nside, x, y)
        counts[pix] = get(counts, pix, 0) + 1
    end
    return counts
end

function _healpix_footprint_weights(nside, lon, lat, n_override)
    spec = HEALPixGridSpec(; nside, T=Float64)
    shape = accumulator_shape(spec)
    # Store only the occupied cells, including at nside=2048.
    data = reshape(spzeros(prod(shape)), shape..., 1)
    weights = spzeros(prod(shape))
    observations = (lat_bounds=reshape(collect(lat), 1, 4),
                    lon_bounds=reshape(collect(lon), 1, 4), values=ones(1, 1))
    SatelliteGridding._accumulate_l2_spherical!(SatelliteGridding._spherical_subdivision(spec),
        data, nothing, reshape(weights, shape), false, spec, observations, n_override)
    return weights
end

function _healpix_hot_loop_allocations(::Type{T}, n_override=nothing) where {T}
    spec = HEALPixGridSpec(nside=8)
    shape = accumulator_shape(spec)
    data, m2, weights = zeros(Float32, shape..., 1), zeros(Float32, shape..., 1), zeros(Float32, shape)
    lon = Matrix{T}(repeat(Int16[177 179 -179 179], 3, 1))
    lat = Matrix{T}(repeat(Int16[-4 -4 4 4], 3, 1))
    capacity = something(n_override, 24)
    values, points = ones(Float32, 3, 1), zeros(capacity, capacity, 2)
    cells = Vector{Tuple{Int,Int}}(undef, capacity^2)
    args = (data, m2, weights, true, spec, lat, lon, values, n_override, points, cells)
    SatelliteGridding._accumulate_unit_vector_footprints!(args...)
    return @allocated SatelliteGridding._accumulate_unit_vector_footprints!(args...)
end

struct _HemisphereTestGrid <: AbstractSphericalGridSpec{Float32} end
SatelliteGridding._spherical_cell(::_HemisphereTestGrid, lat, lon) = (1, lat < 0 ? 1 : 2)
SatelliteGridding.grid_shape(::_HemisphereTestGrid) = (1, 2)
SatelliteGridding._output_spatial_shape(::_HemisphereTestGrid) = (2,)
SatelliteGridding.payload_dims(::_HemisphereTestGrid) = ("hemisphere",)
function SatelliteGridding.write_coordinates!(ds, ::_HemisphereTestGrid)
    defDim(ds, "hemisphere", 2)
    defVar(ds, "hemisphere", Int32, ("hemisphere",))[:] = [-1, 1]
    return nothing
end

@testset "HEALPix grid" begin
    @testset "healpy 1.17.3 NESTED reference" begin
        for (side, longitude, latitude, index) in _healpix_ref_rows("healpix_ang2pix_ref.csv")
            nside, pix = parse(Int, side), parse(Int64, index)
            lon, lat = parse(Float64, longitude), parse(Float64, latitude)
            @test ang2pix_nest(nside, lon, lat) === pix
        end
        for (side, index, longitude, latitude) in _healpix_ref_rows("healpix_pix2ang_ref.csv")
            nside, pix = parse(Int, side), parse(Int64, index)
            lon, lat = pix2ang_nest(nside, pix)
            @test _healpix_angle_error(lon, lat, parse(Float64, longitude), parse(Float64, latitude)) <= 1e-9
            @test ang2pix_nest(nside, lon, lat) == pix
        end
        for (side, index, vertex, longitude, latitude) in _healpix_ref_rows("healpix_corners_ref.csv")
            nside, pix, k = parse(Int, side), parse(Int64, index), parse(Int, vertex)
            lons, lats = healpix_corners(nside, pix)
            @test _healpix_angle_error(lons[k], lats[k], parse(Float64, longitude), parse(Float64, latitude)) <= 1e-9
        end
    end

    @testset "Hierarchy, round trips, and coverage" begin
        rng = MersenneTwister(48)
        for nside in (1, 2, 8, 64, 8192)
            npix = nside2npix(nside)
            @test npix === Int64(12 * nside^2)
            for parent in rand(rng, 0:(npix - 1), 25), child in (4parent):(4parent + 3)
                lon, lat = pix2ang_nest(2nside, child)
                recovered = ang2pix_nest(2nside, lon, lat)
                @test recovered == child
                @test healpix_parent(recovered) == parent
                @test ang2pix_nest(nside, lon, lat) == parent
            end
        end

        nside = 64
        npix = nside2npix(nside)
        @test all(pix -> ang2pix_nest(nside, pix2ang_nest(nside, pix)...) == pix, 0:(npix - 1))
        @test all(0 <= ang2pix_nest(nside, lon, lat) < npix
                  for lat in -90:1.75:90, lon in -720:1.9:720)
        # Edge-line rounding can reach the next sector where the 0° longitude
        # wrap meets the polar/equatorial transition.
        transition = asind(2 / 3)
        for lat in (transition, -transition), adjacent_lat in (prevfloat(lat), lat, nextfloat(lat)),
                lon in (-4e-14, 0.0, 4e-14, prevfloat(360.0))
            pix = ang2pix_nest(nside, lon, adjacent_lat)
            @test 0 <= pix < npix
            centre = pix2ang_nest(nside, pix)
            @test _healpix_angle_error(centre..., lon, adjacent_lat) < 2
        end
        for lat in (-90, -42, 0, 42, 90), lon in (-180, -137.25, 0, 45, 175.75)
            @test ang2pix_nest(nside, lon, lat) == ang2pix_nest(nside, lon - 720, lat)
            @test ang2pix_nest(nside, lon, lat) == ang2pix_nest(nside, lon + 720, lat)
        end
        @test ang2pix_nest(8, 22.5, 0) == ang2pix_nest(8, -337.5, 0) == 278
        # Exercise the full Int64 range supported by the geometry, including
        # cells immediately adjacent to the north and south poles.
        nside = 1 << 29
        npix = nside2npix(nside)
        for pix in (0, nside^2 - 1, 8nside^2, npix - 1, rand(rng, 0:(npix - 1), 12)...)
            @test ang2pix_nest(nside, pix2ang_nest(nside, pix)...) === pix
        end
        spec = HEALPixGridSpec(; nside=2048)
        @test isbits(spec)
        @test accumulator_shape(spec) == (2048^2, 12)
        @test spec.level == 11
        @test SatelliteGridding._storage_kwargs(NetCDFOptions(), spec).chunksizes == [2^20, 1]
        linear_index = LinearIndices(accumulator_shape(spec))
        for pix in rand(rng, 0:(nside2npix(spec.nside) - 1), 12)
            lon, lat = pix2ang_nest(spec.nside, pix)
            col, row = SatelliteGridding._spherical_cell(spec, lat, lon)
            @test linear_index[col, row] == pix + 1
        end
    end

    @testset "Spherical footprint subdivision" begin
        lat = (-4.0, -4.0, 4.0, 4.0)
        lon = (177.0, 179.0, 181.0, 179.0)
        wrapped = (177.0, 179.0, -179.0, 179.0)
        weights = _healpix_footprint_weights(64, lon, lat, 8)
        @test sum(weights) == 1
        @test weights == _healpix_footprint_weights(64, wrapped, lat, 8)
        for order in ((4, 3, 2, 1), (2, 4, 1, 3))
            @test weights == _healpix_footprint_weights(64, map(i -> wrapped[i], order),
                                                        map(i -> lat[i], order), 8)
        end

        lon = (44.487, 45.349, 31.818, 31.034)
        lat = (87.116, 87.171, 87.829, 87.774)
        @test all(ang2pix_nest(64, x, y) == 4092 for (x, y) in zip(lon, lat))
        # The review's 0.75/0.25 split comes from lon/lat subdivision.
        # The geodesic footprint of the same corners stays in cell 4092.
        planar = _healpix_reference_counts(64, lon, lat, 8; spherical=false)
        @test planar == Dict(4092 => 48, 4089 => 16)
        dense_planar = _healpix_reference_counts(64, lon, lat, 512; spherical=false)
        @test dense_planar[4089] / 512^2 ≈ 0.25 atol=0.02
        dense_spherical = _healpix_reference_counts(64, lon, lat, 512)
        @test dense_spherical == Dict(4092 => 512^2)
        @test findnz(_healpix_footprint_weights(64, lon, lat, 8)) == ([4093], [1.0])

        # This footprint crosses a curved edge even in 3-D, despite all four
        # corners being in one cell. A corner shortcut would lose cell 4034.
        lon = (41.994000000000014, 38.57234680759891, 38.489880273660205, 41.910179640718574)
        lat = (79.01871394891049, 79.75136013762997, 79.7293429788161, 78.99668369523602)
        @test all(ang2pix_nest(64, x, y) == 3991 for (x, y) in zip(lon, lat))
        dense = _healpix_reference_counts(64, lon, lat, 512)
        cells, fractions = findnz(_healpix_footprint_weights(64, lon, lat, 8))
        @test cells .- 1 == sort!(collect(keys(dense))) == [3991, 4034]
        @test fractions == [61, 3] ./ 64
        @test fractions[2] ≈ dense[4034] / 512^2 atol=1/64

        lon, lat = (-135.0, -45.0, 45.0, 135.0), (89.0, 89.0, 89.0, 89.0)
        nside = 128
        polar = _healpix_footprint_weights(nside, lon, lat, 8)
        @test [sum(polar[face*nside^2 .+ (1:nside^2)]) for face in 0:3] == fill(0.25, 4)
        @test sum(polar) == 1

        # An odd number of samples lands exactly on the enclosed pole. Its
        # longitude is defined as 0° and all other samples split symmetrically.
        points = zeros(9, 9, 2)
        corners, normal = SatelliteGridding._sphere_ordered_corners(
            map(SatelliteGridding._sphere_unit_vector, lat, lon)...)
        SatelliteGridding._sphere_subpoints!(points, corners, normal, 9)
        @test points[5, 5, :] == [90.0, 0.0]
        polar = _healpix_footprint_weights(nside, lon, lat, 9)
        polar[ang2pix_nest(nside, 0, 90) + 1] -= 1/81
        @test [sum(polar[face*nside^2 .+ (1:nside^2)]) for face in 0:3] ≈ fill(20/81, 4)
        @test _healpix_hot_loop_allocations(Float64) == 0
        @test _healpix_hot_loop_allocations(Union{Int16,Float64}) == 0
        @test _healpix_hot_loop_allocations(Float64, 31) == 0
    end

    @testset "Single-cell HEALPix footprints retain integral weights" begin
        for (T, n_override, n_footprints) in ((Float32, 256, 257), (Float16, nothing, 129))
            spec = HEALPixGridSpec(nside=8; T)
            lon, lat = pix2ang_nest(spec.nside, 0)
            lon_bounds = repeat([lon - 0.01 lon + 0.01 lon + 0.01 lon - 0.01], n_footprints, 1)
            lat_bounds = repeat([lat - 0.01 lat - 0.01 lat + 0.01 lat + 0.01], n_footprints, 1)
            observations = (; lat_bounds, lon_bounds, values=fill(T(3), n_footprints, 1))
            shape = accumulator_shape(spec)
            data, m2, weights = zeros(T, shape..., 1), zeros(T, shape..., 1), zeros(T, shape)
            SatelliteGridding._accumulate_l2_spherical!(
                SatelliteGridding._spherical_subdivision(spec), data, m2, weights,
                true, spec, observations, n_override)
            @test sum(weights) === T(n_footprints)
            @test weights[1, 1] === T(n_footprints)
            @test data[1, 1, 1] === T(3)
            @test all(iszero, m2)
        end
    end

    @testset "Adaptive footprint oversampling and fixed override" begin
        # 14.5 × 5.5 km near the equator. Dense sampling hits 14 cells, four
        # of which would have been missed by the former 4 × 4 default.
        nside = 2048
        width, height = rad2deg(14_500 / 6.371e6), rad2deg(5_500 / 6.371e6)
        lon = (0.001, 0.001 + width, 0.001 + width, 0.001)
        lat = (0.001, 0.001, 0.001 + height, 0.001 + height)
        dense = _healpix_reference_counts(nside, lon, lat, 512)
        adaptive = _healpix_footprint_weights(nside, lon, lat, nothing)
        fixed = _healpix_footprint_weights(nside, lon, lat, 4)
        @test findnz(adaptive)[1] .- 1 == sort!(collect(keys(dense)))
        @test adaptive == _healpix_footprint_weights(nside, lon, lat, 20)
        @test length(dense) == 14
        @test nnz(fixed) == 10
        @test sum(adaptive) ≈ 1

        width, height = rad2deg(2_250 / 6.371e6), rad2deg(1_300 / 6.371e6)
        lon = (-width/2, width/2, width/2, -width/2)
        lat = (-height/2, -height/2, height/2, height/2)
        spec = HEALPixGridSpec(nside=64)
        corners = map(_healpix_unit_vector, lon, lat)
        @test SatelliteGridding._unit_vector_samples(corners,
            SatelliteGridding._healpix_cell_spacing(spec)) == 4
        @test _healpix_footprint_weights(64, lon, lat, nothing) ==
              _healpix_footprint_weights(64, lon, lat, 4)
    end

    @testset "CLI and invalid inputs" begin
        args = parse_l2_args(["--config", "input.toml", "--gridType", "healpix", "--nside", "64"])
        spec = args_to_grid_spec(args)
        @test spec isa HEALPixGridSpec{Float32}
        @test spec.nside == 64
        @test spec.level == 6
        @test !spec.cell_bounds
        @test resolve_backend("sequential"; grid_spec=spec) === nothing
        @test_throws "sequential CPU" resolve_backend("cuda"; grid_spec=spec)
        @test_throws "sequential CPU" resolve_backend("metal"; grid_spec=spec)
        for invalid in (0, -8, 3, 96, 1 << 30)
            args["nside"] = invalid
            @test_throws ArgumentError args_to_grid_spec(args)
        end
        @test_throws ArgumentError HEALPixGridSpec(nside=4.0)
        @test_throws ArgumentError HEALPixGridSpec(nside=true)
        @test_throws ArgumentError nside2npix(0)
        @test_throws ArgumentError ang2pix_nest(2, Inf, 0)
        @test_throws ArgumentError ang2pix_nest(2, 0, -91)
        @test_throws ArgumentError pix2ang_nest(2, -1)
        @test_throws ArgumentError healpix_corners(2, nside2npix(2))
    end

    @testset "Oversampling validation before output creation" begin
        mktempdir() do dir
            config = DataSourceConfig(Dict{String,String}(), OrderedDict("signal" => "signal"),
                                      FilterRule[], "missing_YYYY-MM-DD.nc", dir)
            date = DateTime(2020, 1, 1)
            times = TimeSpec(date, date, Day(1))
            output = joinpath(dir, "invalid_oversampling.nc")
            @test SubpixelGridding(Int32(7)).n_oversample === 7
            for invalid in (0, -1, 4.0, 4.5, true, "4")
                @test_throws ArgumentError SubpixelGridding(invalid)
                @test_throws ArgumentError SubpixelGridding(n_oversample=invalid)
                @test_throws ArgumentError grid_l2(config, HEALPixGridSpec(nside=8), times;
                                                   n_oversample=invalid, outfile=output)
                @test_throws ArgumentError grid_l2(config, HEALPixGridSpec(nside=8), times;
                    n_oversample=invalid, footprint_method=SubpixelGridding(8), outfile=output)
                @test !ispath(output)
            end
        end
    end

    @testset "Packed centres and bounds give the same weights as physical degrees" begin
        mktempdir() do dir
            basic = Dict("lat" => "geo/lat", "lon" => "geo/lon",
                         "lat_bnd" => "geo/lat_bounds", "lon_bnd" => "geo/lon_bounds")
            config = DataSourceConfig(basic, OrderedDict("signal" => "signal"),
                                      FilterRule[], "packed_YYYY-MM-DD.nc", dir)
            input = joinpath(dir, "packed_2020-01-01.nc")
            output = joinpath(dir, "packed_output.nc")
            date = DateTime(2020, 1, 1)
            times = TimeSpec(date, date, Day(1))
            grids = (RectangularGridSpec(dlat=2f0, dlon=2f0),
                     CubedSphereGridSpec(Nc=8), HEALPixGridSpec(nside=8))
            centres = (Int16[150, 250, 350], Int16[550, 1550, 3550])
            corners = (Int16[100 100 200 200; 200 200 300 300; 300 300 400 400],
                       Int16[500 600 600 500; 1500 1600 1600 1500; 3500 3600 3600 3500])
            expected = Dict()
            for packing in (:unpacked, :centres, :bounds, :both)
                Dataset(input, "c") do ds
                    defDim(ds, "sounding", 3)
                    defDim(ds, "vertex", 4)
                    geo = defGroup(ds, "geo")
                    for (names, raw, dims, packed) in (
                        (("lat", "lon"), centres, ("sounding",), packing in (:centres, :both)),
                        (("lat_bounds", "lon_bounds"), corners, ("sounding", "vertex"), packing in (:bounds, :both)),
                    )
                        for (name, values, offset) in zip(names, raw, (0.0, 4.0))
                            selectors = ntuple(_ -> Colon(), length(dims))
                            if packed
                                var = defVar(geo, name, Int16, dims;
                                    attrib=["scale_factor" => 0.01, "add_offset" => offset])
                                var.var[selectors...] = values
                            else
                                var = defVar(geo, name, Float64, dims)
                                var.var[selectors...] = values .* 0.01 .+ offset
                            end
                        end
                    end
                    defVar(ds, "signal", Float32, ("sounding",))[:] = [10, 20, 30]
                end
                for spec in grids
                    grid_l2(config, spec, times; n_oversample=8, outfile=output, progress=false)
                    Dataset(output) do ds
                        weights = ds["n"][:]
                        @test sum(weights) == 3
                        if packing === :unpacked
                            expected[typeof(spec)] = weights
                        else
                            @test weights == expected[typeof(spec)]
                        end
                    end
                end
            end
        end
    end

    @testset "Invalid geolocation preserves the rest of the input file" begin
        mktempdir() do dir
            input = joinpath(dir, "geo_2020-01-01.nc")
            output = joinpath(dir, "geo_output.nc")
            config = DataSourceConfig(Dict("lat" => "lat", "lon" => "lon",
                                          "lat_bnd" => "lat_bounds", "lon_bnd" => "lon_bounds"),
                                      OrderedDict("signal" => "signal"), FilterRule[],
                                      "geo_YYYY-MM-DD.nc", dir)
            date = DateTime(2020, 1, 1)
            times = TimeSpec(date, date, Day(1))
            cases = ((:fill, -9999f0), (:fill, 2.5f0), (:integer_fill, Int16(-9999)),
                     (:raw_lat, -9999f0), (:missing_lon, 15f0), (:missing_lon, -0f0),
                     (:centre_lat, 91f0), (:centre_lon, Inf32), (:corner_lon, Inf32),
                     (:centre_lon, 720.01f0), (:corner_lon, -9999f0),
                     (:default_centre_lon, 9.969209968386869e36), (:default_corner_lon, 9.96921f36),
                     (:default_packed_lon, Int16(-32767)),
                     (:packed_fill, Int16(-800)), (:packed_missing, Int16(900)),
                     (:polar_overshoot, nextfloat(90f0)),
                     (:polar_overshoot, Float64(nextfloat(90f0))),
                     (:polar_overshoot, -Float64(nextfloat(90f0))))
            for (kind, invalid) in cases
                Dataset(input, "c") do ds
                    defDim(ds, "sounding", 3)
                    defDim(ds, "vertex", 4)
                    centre_dims = ("sounding",)
                    corner_dims = ("sounding", "vertex")
                    lat_type = kind === :polar_overshoot ? typeof(invalid) : Float32
                    lat = defVar(ds, "lat", lat_type, centre_dims)
                    lon = defVar(ds, "lon", kind === :default_centre_lon ? Float64 : Float32, centre_dims)
                    lat_bounds_type = kind === :integer_fill ? Int16 : lat_type
                    lat_storage = kind in (:fill, :integer_fill) ? (; fillvalue=invalid) : (;)
                    lat_bounds = defVar(ds, "lat_bounds", lat_bounds_type, corner_dims; lat_storage...)
                    packed_lon = kind in (:default_packed_lon, :packed_fill, :packed_missing)
                    lon_attributes = Dict{String,Any}()
                    packed_lon && (lon_attributes["scale_factor"] = 0.01)
                    kind === :missing_lon && (lon_attributes["missing_value"] = [invalid, -9999f0])
                    kind === :packed_missing && (lon_attributes["missing_value"] = [invalid, Int16(-900)])
                    lon_storage = kind === :packed_fill ? (; fillvalue=invalid) : (;)
                    lon_bounds = defVar(ds, "lon_bounds", packed_lon ? Int16 : Float32, corner_dims;
                        attrib=lon_attributes, lon_storage...)
                    lat[:] = Float32[1.5, 2.5, 3.5]
                    lon.var[1] = 5.5
                    kind === :default_centre_lon || (lon.var[2] = 15.5)
                    lon.var[3] = 35.5
                    lat_bounds.var[:, :] = Float32[1 1 2 2; 2 2 3 3; 3 3 4 4]
                    lon_corners = [5 6 6 5; 15 16 16 15; 35 36 36 35]
                    packed_lon && (lon_corners .*= 100)
                    lon_bounds.var[:, 2:4] = lon_corners[:, 2:4]
                    lon_bounds.var[1, 1] = lon_corners[1, 1]
                    kind in (:default_corner_lon, :default_packed_lon) ||
                        (lon_bounds.var[2, 1] = lon_corners[2, 1])
                    lon_bounds.var[3, 1] = lon_corners[3, 1]
                    defVar(ds, "signal", Float32, centre_dims)[:] = Float32[10, 20, 30]
                    if kind === :polar_overshoot
                        lat[2] = invalid
                        near_pole = copysign(89, invalid)
                        lat_bounds.var[2, :] = [near_pole, near_pole, invalid, invalid]
                    elseif kind === :centre_lat
                        lat[2] = invalid
                    elseif kind === :centre_lon
                        lon[2] = invalid
                    elseif kind === :missing_lon
                        lon_bounds.var[2, 1] = abs(invalid)
                    elseif kind in (:corner_lon, :packed_fill, :packed_missing)
                        lon_bounds.var[2, 1] = invalid
                    elseif kind in (:default_centre_lon, :default_corner_lon, :default_packed_lon)
                        var = kind === :default_centre_lon ? lon : lon_bounds
                        @test !haskey(var.attrib, "_FillValue")
                        @test var.var[2, 1] === invalid
                    else
                        lat_bounds.var[2, 1] = invalid
                    end
                end
                if kind in (:default_centre_lon, :default_corner_lon, :default_packed_lon)
                    Dataset(input) do ds
                        bounds = kind !== :default_centre_lon
                        path = bounds ? "lon_bounds" : "lon"
                        @test isnan(SatelliteGridding._read_l2_geolocation(ds, path; bounds)[2, 1])
                    end
                end
                grids = (HEALPixGridSpec(nside=8), CubedSphereGridSpec(Nc=8))
                # These cases also exercise the rectangular geolocation path.
                if kind in (:fill, :integer_fill, :raw_lat, :missing_lon)
                    grids = (grids..., RectangularGridSpec(dlat=10f0, dlon=10f0))
                end
                for spec in grids
                    grid_l2(config, spec, times; n_oversample=8, outfile=output, progress=false)
                    Dataset(output) do ds
                        expected = kind === :polar_overshoot ? [10f0, 20f0, 30f0] : [10f0, 30f0]
                        @test sum(ds["n"][:]) == length(expected)
                        @test sort!(unique!(collect(skipmissing(ds["signal"][:])))) == expected
                        conventions = spec isa HEALPixGridSpec ? "CF-1.13" :
                                      spec isa CubedSphereGridSpec ? "CF" : nothing
                        @test get(ds.attrib, "Conventions", nothing) == conventions
                    end
                end
            end
        end
    end

    @testset "Minimal spherical grid can write NetCDF" begin
        mktempdir() do dir
            _accuracy_input(joinpath(dir, "demo_2020-01-01.nc"), [(15.0, 16.0, -2.0, 2.0)], Float32[10])
            config = DataSourceConfig(Dict("lat" => "lat", "lon" => "lon",
                                          "lat_bnd" => "lat_bounds", "lon_bnd" => "lon_bounds"),
                                      OrderedDict("signal" => "signal"), FilterRule[],
                                      "demo_YYYY-MM-DD.nc", dir)
            date = DateTime(2020, 1, 1)
            times = TimeSpec(date, date, Day(1))
            output = joinpath(dir, "demo.nc")
            grid(config, _HemisphereTestGrid(), times; outfile=output, progress=false)
            Dataset(output) do ds
                @test ds.attrib["title"] == "Gridded satellite data"
                @test dimnames(ds["signal"]) == ("hemisphere", "time")
                @test !haskey(ds["signal"].attrib, "grid_mapping")
                @test !haskey(ds["signal"].attrib, "coordinates")
                @test ds["hemisphere"][:] == [-1, 1]
                @test ds["signal"][:, 1] == [10, 10]
                @test ds["n"][:, 1] == [0.5, 0.5]
            end
        end
    end

    @testset "L2 footprints to NetCDF" begin
        mktempdir() do dir
            input = joinpath(dir, "HEALPIX_2020-01-01.nc")
            # The base-face junction at (45°, 0°), a dateline crossing, a
            # northern polar-cap footprint, and a footprint rejected by QA.
            base_edge = (40.0, 50.0, -2.0, 2.0)
            _accuracy_input(input, [base_edge, (177.0, -177.0, -3.0, 3.0),
                                    (20.0, 80.0, 86.0, 89.5), base_edge],
                            Float32[10, 20, 30, 99]; qa=Float32[1, 1, 1, 0])
            Dataset(input, "a") do ds
                ds["lon"][2] = 180  # the fixture's arithmetic mean crosses the seam
                ds["signal"].attrib["coordinates"] = "source_lat source_lon"
                ds["signal"].attrib["grid_mapping"] = "source_mapping"
            end
            cfgfile = joinpath(dir, "healpix.toml")
            write(cfgfile, """
                filePattern = "HEALPIX_YYYY-MM-DD.nc"
                folder = "$(replace(dir, "\\" => "/"))/"

                [basic]
                lat = "lat"
                lon = "lon"
                lat_bnd = "lat_bounds"
                lon_bnd = "lon_bounds"

                [grid]
                signal = "signal"

                [filter]
                qa = "> 0"

                [metadata]
                variable_attributes = ["units", "coordinates", "grid_mapping"]
                """)
            config = load_config(cfgfile)
            date = DateTime("2020-01-01")
            times = TimeSpec(date, date, Day(1))
            spec = HEALPixGridSpec(nside=8, cell_bounds=true)
            npix = nside2npix(spec.nside)
            output = joinpath(dir, "healpix.nc")
            grid(config, spec, times, SubpixelGridding(n_oversample=8);
                 outfile=output, compute_std=true, progress=false)

            Dataset(output) do ds
                @test ds.dim["cell"] == npix
                @test ds.dim["time"] == 1
                @test ds.dim["nv"] == 4
                @test ds.attrib["Conventions"] == "CF-1.13"
                @test ds.attrib["grid"] == "healpix"
                @test !haskey(ds.attrib, "grid_mapping_name")
                @test ds.attrib["created_with"] == "SatelliteGridding.jl"
                @test ds.attrib["nside"] == 8
                @test ds.attrib["healpix_level"] == 3
                @test ds.attrib["healpix_ordering"] == "nested"
                ids = ds["cell_ids"]
                @test dimnames(ids) == ("cell",)
                @test eltype(ids.var) == Int64
                @test ids[:] == 0:(npix - 1)
                @test ids.attrib["standard_name"] == "healpix_index"
                @test ids.attrib["units"] == "1"
                @test ids.attrib["grid_name"] == "healpix"
                @test ids.attrib["ellipsoid"] == "sphere"
                @test ids.attrib["level"] == 3
                @test ids.attrib["indexing_scheme"] == "nested"
                mapping = ds["healpix"]
                @test dimnames(mapping) == ()
                @test eltype(mapping.var) == Int32
                @test mapping.attrib["grid_mapping_name"] == "healpix"
                @test mapping.attrib["refinement_level"] == 3
                @test mapping.attrib["indexing_scheme"] == "nested"
                @test mapping.attrib["earth_radius"] == 6.371e6
                area = ds["cell_area"]
                @test dimnames(area) == ()
                @test area[] ≈ 4pi * (6.371e6)^2 / npix
                @test area.attrib["units"] == "m2"
                @test occursin("equal solid angle", area.attrib["comment"])
                @test occursin("1.35%", area.attrib["comment"])

                for (name, standard_name, units) in (("lon", "longitude", "degrees_east"),
                                                       ("lat", "latitude", "degrees_north"))
                    coord = ds[name]
                    @test dimnames(coord) == ("cell",)
                    @test eltype(coord.var) == Float64
                    @test coord.attrib["standard_name"] == standard_name
                    @test coord.attrib["units"] == units
                    @test coord.attrib["bounds"] == name * "_bounds"
                    @test dimnames(ds[name * "_bounds"]) == ("nv", "cell")
                    @test eltype(ds[name * "_bounds"].var) == Float64
                    @test occursin("Corner vertices for plotting", ds[name * "_bounds"].attrib["comment"])
                end
                @test (ds["lon"][end], ds["lat"][end]) == pix2ang_nest(spec.nside, npix - 1)
                corners = healpix_corners(spec.nside, npix - 1)
                @test Tuple(ds["lon_bounds"][:, end]) == corners[1]
                @test Tuple(ds["lat_bounds"][:, end]) == corners[2]

                for name in ("signal", "signal_std", "n")
                    var = ds[name]
                    @test dimnames(var) == ("cell", "time")
                    @test size(var) == (npix, 1)
                    @test var.attrib["coordinates"] == "lat lon cell_ids"
                    @test var.attrib["grid_mapping"] == "healpix"
                    @test Tuple(NCDatasets.chunking(var.var)[2]) == (npix, 1)
                end
                weights, values = ds["n"][:, 1], ds["signal"][:, 1]
                @test ds["signal"].attrib["units"] == "ppm"
                @test ds["n"].attrib["long_name"] == "Sum of fractional footprint overlap weights"
                @test sum(weights) == 3
                @test ismissing.(values) == iszero.(weights)
                @test all(iszero, skipmissing(ds["signal_std"][:, 1]))
                hits(value) = findall(v -> isequal(v, value), values)
                base, seam, polar = hits(10f0), hits(20f0), hits(30f0)
                @test length(unique((base .- 1) .÷ spec.nside^2)) >= 2
                @test length(seam) >= 2
                seam_lons = ds["lon"][:][seam]
                @test minimum(seam_lons) < 180 < maximum(seam_lons)
                @test all(abs.(seam_lons .- 180) .< 15)
                @test !isempty(polar) && all(ds["lat"][:][polar] .> 80)
                @test sum(weights[base]) == sum(weights[seam]) == sum(weights[polar]) == 1
            end

            rolling = TimeSpec(date, date, Day(1); sample_step=Day(1), window_halfwidth_days=1)
            rolled = joinpath(dir, "healpix_rolling.nc")
            grid_l2(config, HEALPixGridSpec(nside=8, T=Float64), rolling;
                    outfile=rolled, n_oversample=8, progress=false)
            Dataset(rolled) do ds
                @test dimnames(ds["n"]) == ("cell", "time")
                @test sum(ds["n"][:, 1]) == 3
                @test ds["window_start"][1] == date - Day(1)
                @test ds["window_end"][1] == date + Day(1)
                @test !haskey(ds.dim, "nv")
                @test !haskey(ds, "lon_bounds") && !haskey(ds, "lat_bounds")
                @test !haskey(ds["lon"].attrib, "bounds")
            end

            @test_throws "sequential CPU" grid_l2(config, spec, times;
                                                   backend=SatelliteGridding.CPU(), outfile=output)
            @test_throws "CircularFootprintGridding" grid(config, spec, times, CircularFootprintGridding())
            @test_throws "rectangular" grid_center(config, spec, times)
            @test_throws "rectangular" grid(config, spec, times, CenterPointGridding())
            @test_throws "Zarr output requires" grid_l2(config, spec, times; zarr_out=joinpath(dir, "x.zarr"))
            @test_throws "Zarr output requires" init_zarr_store(joinpath(dir, "x.zarr"), spec, config.grid_vars)
            @test_throws "Zarr output requires" write_zarr_meta_json(dir, spec, config.grid_vars)
        end
    end

    @testset "Coordinate writing across chunks" begin
        mktempdir() do dir
            spec = HEALPixGridSpec(nside=128, cell_bounds=true)
            path = joinpath(dir, "chunks.nc")
            ds, vars = create_output_dataset(path, spec, 1, OrderedDict("x" => "x"), false;
                output_options=NetCDFOptions(chunksizes=(100, 1)))
            @test Tuple(NCDatasets.chunking(vars["x"].var)[2]) == (100, 1)
            close(ds)
            Dataset(path) do ds
                chunks = NCDatasets.chunking(ds["lon"].var)[2]
                @test only(chunks) < ds.dim["cell"]
                for cell in (only(chunks), only(chunks) + 1, ds.dim["cell"])
                    @test ds["cell_ids"][cell] == cell - 1
                    @test (ds["lon"][cell], ds["lat"][cell]) == pix2ang_nest(spec.nside, cell - 1)
                    lons, lats = healpix_corners(spec.nside, cell - 1)
                    @test Tuple(ds["lon_bounds"][:, cell]) == lons
                    @test Tuple(ds["lat_bounds"][:, cell]) == lats
                end
            end
        end
    end
end
