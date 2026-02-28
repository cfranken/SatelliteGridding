@testset "Vegetation Indices" begin

    @testset "NDVI" begin
        # red=0.1, nir=0.5 → NDVI = (0.5-0.1)/(0.5+0.1) = 0.4/0.6 ≈ 0.6667
        @test compute_ndvi(0.1f0, 0.5f0) ≈ 0.6667f0 rtol=1e-3

        # Equal red and nir → NDVI = 0
        @test compute_ndvi(0.3f0, 0.3f0) ≈ 0.0f0

        # NIR >> RED → NDVI ≈ 1
        @test compute_ndvi(0.01f0, 1.0f0) > 0.95f0
    end

    @testset "EVI" begin
        red = 0.1f0
        nir = 0.5f0
        blue = 0.05f0
        evi = compute_evi(red, nir, blue)

        # EVI = 2.5 * (0.5-0.1) / (0.5 + 6*0.1 - 7.5*0.05 + 10000)
        expected = 2.5f0 * 0.4f0 / (0.5f0 + 0.6f0 - 0.375f0 + 10000.0f0)
        @test evi ≈ expected rtol=1e-4
    end

    @testset "NIRv" begin
        red = 0.1f0
        nir = 0.5f0
        nirv = compute_nirv(red, nir)
        expected_ndvi = (0.5f0 - 0.1f0) / (0.5f0 + 0.1f0)
        @test nirv ≈ expected_ndvi * nir rtol=1e-5
    end

    @testset "NDWI" begin
        nir = 0.5f0
        swir = 0.2f0
        ndwi = compute_ndwi(nir, swir)
        @test ndwi ≈ (0.5f0 - 0.2f0) / (0.5f0 + 0.2f0) rtol=1e-5
    end
end
