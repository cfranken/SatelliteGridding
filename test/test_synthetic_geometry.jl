function _exact_rectangle_weights(xmin, xmax, ymin, ymax; n_lon=6, n_lat=6)
    lat_corners = Float32[ymin, ymin, ymax, ymax]
    lon_corners = Float32[xmin, xmax, xmax, xmin]
    footprint_weights(ExactIntersectionGridding(), lat_corners, lon_corners;
                      n_lon=n_lon, n_lat=n_lat)
end

function _subpixel_rectangle_weights(xmin, xmax, ymin, ymax, n_oversample;
                                     n_lon=6, n_lat=6)
    lat_corners = Float32[ymin, ymin, ymax, ymax]
    lon_corners = Float32[xmin, xmax, xmax, xmin]
    footprint_weights(SubpixelGridding(n_oversample), lat_corners, lon_corners;
                      n_lon=n_lon, n_lat=n_lat)
end

function _accumulate_synthetic_rectangles(rectangles, values, n_oversample; n_lon=6, n_lat=6)
    n_pixels = length(rectangles)
    grid_data = zeros(Float32, n_lon, n_lat, 1)
    grid_std = zeros(Float32, n_lon, n_lat, 1)
    grid_weights = zeros(Float32, n_lon, n_lat)

    points = zeros(Float32, n_oversample, n_oversample, 2)
    ix = zeros(Int32, n_oversample^2)
    iy = zeros(Int32, n_oversample^2)
    lats0 = zeros(Float32, n_oversample)
    lons0 = zeros(Float32, n_oversample)
    lats1 = zeros(Float32, n_oversample)
    lons1 = zeros(Float32, n_oversample)

    lat_idx = zeros(Float32, n_pixels, 4)
    lon_idx = zeros(Float32, n_pixels, 4)
    for (i, (xmin, xmax, ymin, ymax)) in enumerate(rectangles)
        lat_idx[i, :] = Float32[ymin, ymin, ymax, ymax]
        lon_idx[i, :] = Float32[xmin, xmax, xmax, xmin]
    end

    accumulate_footprint!(grid_data, grid_std, grid_weights, false,
                          lat_idx, lon_idx, reshape(Float32.(values), :, 1),
                          n_pixels, 1, n_oversample,
                          points, ix, iy, lats0, lons0, lats1, lons1)
    grid_data[:, :, 1], grid_weights
end

function _exact_weighted_mean(rectangles, values; n_lon=6, n_lat=6)
    weighted_sum = zeros(Float32, n_lon, n_lat)
    weights = zeros(Float32, n_lon, n_lat)
    for (rect, value) in zip(rectangles, values)
        w = _exact_rectangle_weights(rect...; n_lon=n_lon, n_lat=n_lat)
        weights .+= w
        weighted_sum .+= Float32(value) .* w
    end

    mean = zeros(Float32, n_lon, n_lat)
    for i in eachindex(weights)
        weights[i] > 0 && (mean[i] = weighted_sum[i] / weights[i])
    end
    mean, weights
end

@testset "Synthetic Geometry Baselines" begin
    @testset "Subpixel rectangle weights converge toward exact area overlap" begin
        rect = (1.17, 4.73, 1.29, 3.61)
        expected = _exact_rectangle_weights(rect...)

        w5 = _subpixel_rectangle_weights(rect..., 5)
        w20 = _subpixel_rectangle_weights(rect..., 20)
        w100 = _subpixel_rectangle_weights(rect..., 100)

        err5 = maximum(abs.(w5 .- expected))
        err20 = maximum(abs.(w20 .- expected))
        err100 = maximum(abs.(w100 .- expected))

        @test sum(w5) ≈ 1.0f0 rtol=1f-5
        @test sum(w100) ≈ 1.0f0 rtol=2f-5
        @test err20 < err5
        @test err100 < err20
        @test err100 < 0.003f0
    end

    @testset "Weighted means match exact overlap for simple fake footprints" begin
        rectangles = [
            (1.17, 4.73, 1.29, 3.61),
            (2.05, 4.45, 2.10, 4.20),
        ]
        values = Float32[10, 30]

        expected_mean, expected_weights = _exact_weighted_mean(rectangles, values)
        actual_mean, actual_weights = _accumulate_synthetic_rectangles(rectangles, values, 100)

        mask = expected_weights .> 0.02f0
        @test maximum(abs.(actual_weights .- expected_weights)) < 0.008f0
        @test maximum(abs.(actual_mean[mask] .- expected_mean[mask])) < 0.15f0
    end

    @testset "Method dispatch exposes approximate, center, and exact hooks" begin
        @test SubpixelGridding().n_oversample === nothing
        @test SubpixelGridding(12).n_oversample == 12
        @test CenterPointGridding() isa AbstractGriddingMethod
        @test ExactIntersectionGridding() isa AbstractGriddingMethod

        config = DataSourceConfig(Dict{String,String}(),
                                  OrderedDict("value" => "value"),
                                  FilterRule[], "*.nc", "/tmp")
        grid_spec = GridSpec(lat_min=0.0f0, lat_max=1.0f0,
                             lon_min=0.0f0, lon_max=1.0f0,
                             dlat=1.0f0, dlon=1.0f0)
        time_spec = TimeSpec(DateTime("2020-01-01"), DateTime("2020-01-01"),
                             Dates.Day(1))
        @test_throws ErrorException grid(config, grid_spec, time_spec,
                                         ExactIntersectionGridding())

        lat_corners = Float32[1.25, 1.25, 2.75, 2.75]
        lon_corners = Float32[1.25, 2.75, 2.75, 1.25]
        exact = footprint_weights(ExactIntersectionGridding(), lat_corners, lon_corners;
                                  n_lon=4, n_lat=4)
        center = footprint_weights(CenterPointGridding(), lat_corners, lon_corners;
                                   n_lon=4, n_lat=4)
        @test sum(exact) ≈ 1.0f0
        @test center[2, 2] == 1.0f0
        @test sum(center) == 1.0f0
    end
end
