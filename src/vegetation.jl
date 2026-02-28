"""
    compute_evi(red, nir, blue; G=2.5f0, C1=6.0f0, C2=7.5f0, L=10000.0f0)

Compute the Enhanced Vegetation Index (EVI).

    EVI = G × (NIR - RED) / (NIR + C1×RED - C2×BLUE + L)
"""
function compute_evi(red::T, nir::T, blue::T;
                     G::T=T(2.5), C1::T=T(6.0), C2::T=T(7.5), L::T=T(10000.0)) where {T<:Real}
    G * (nir - red) / (nir + C1 * red - C2 * blue + L)
end

"""
    compute_ndvi(red, nir)

Compute the Normalized Difference Vegetation Index (NDVI).

    NDVI = (NIR - RED) / (NIR + RED)
"""
function compute_ndvi(red::T, nir::T) where {T<:Real}
    (nir - red) / (nir + red)
end

"""
    compute_nirv(red, nir)

Compute the Near-Infrared Reflectance of Vegetation (NIRv).

    NIRv = NDVI × NIR
"""
function compute_nirv(red::T, nir::T) where {T<:Real}
    compute_ndvi(red, nir) * nir
end

"""
    compute_ndwi(nir, swir)

Compute the Normalized Difference Water Index (NDWI).

    NDWI = (NIR - SWIR) / (NIR + SWIR)
"""
function compute_ndwi(nir::T, swir::T) where {T<:Real}
    (nir - swir) / (nir + swir)
end
