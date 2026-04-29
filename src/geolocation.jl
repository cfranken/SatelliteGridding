const MODIS_SINUSOIDAL_RADIUS = 6371007.181
const MODIS_HORIZONTAL_TILES = 36
const MODIS_VERTICAL_TILES = 18
const MODIS_DEFAULT_PIXELS = 2400

abstract type AbstractCenterGeolocation end

struct BasicVariableGeolocation <: AbstractCenterGeolocation end

mutable struct ModisMonolithicLUTGeolocation <: AbstractCenterGeolocation
    path::String
    dataset::Any
end

ModisMonolithicLUTGeolocation(path::String) = ModisMonolithicLUTGeolocation(path, nothing)

struct ModisSinusoidalGeolocation <: AbstractCenterGeolocation
    cache_dir::String
    pixels::Int
end

"""
    default_cache_dir() -> String

Return the base cache directory for generated auxiliary data. Override with
`SATELLITEGRIDDING_CACHE`.
"""
default_cache_dir() =
    get(ENV, "SATELLITEGRIDDING_CACHE",
        joinpath(homedir(), ".cache", "SatelliteGridding"))

"""
    default_modis_cache_dir(; pixels=2400) -> String

Return the default per-tile MODIS sinusoidal geolocation cache directory.
Override with `SATELLITEGRIDDING_MODIS_GEO_CACHE`.
"""
function default_modis_cache_dir(; pixels::Int=MODIS_DEFAULT_PIXELS)
    default_dir = joinpath(default_cache_dir(), "modis", "sinusoidal_$(pixels)px_v1")
    get(ENV, "SATELLITEGRIDDING_MODIS_GEO_CACHE", default_dir)
end

function _validate_modis_tile(h::Integer, v::Integer)
    0 <= h < MODIS_HORIZONTAL_TILES ||
        error("MODIS horizontal tile h must be in 0:35, got $h")
    0 <= v < MODIS_VERTICAL_TILES ||
        error("MODIS vertical tile v must be in 0:17, got $v")
    nothing
end

function _validate_modis_pixels(pixels::Integer)
    pixels > 0 || error("MODIS pixels per tile must be positive, got $pixels")
    nothing
end

"""
    parse_modis_tile(path_or_tile) -> (h, v)

Parse a MODIS tile identifier such as `h08v04` from a filename or string.
"""
function parse_modis_tile(path_or_tile::AbstractString)
    m = match(r"(?:^|[.])h(\d{2})v(\d{2})(?:[.]|$)", basename(String(path_or_tile)))
    m === nothing && (m = match(r"^h(\d{2})v(\d{2})$", String(path_or_tile)))
    m === nothing && error("Could not parse MODIS h/v tile from: $path_or_tile")
    h = parse(Int, m.captures[1])
    v = parse(Int, m.captures[2])
    _validate_modis_tile(h, v)
    h, v
end

modis_tile_width() = 2 * pi * MODIS_SINUSOIDAL_RADIUS / MODIS_HORIZONTAL_TILES
modis_pixel_size(pixels::Integer=MODIS_DEFAULT_PIXELS) = modis_tile_width() / pixels

"""
    modis_sinusoidal_latlon(h, v, line, sample; pixels=2400)

Compute the latitude/longitude at a MODLAND sinusoidal pixel center. `h` and
`v` are zero-based tile indices; `line` and `sample` are one-based pixel
indices within the tile.
"""
function modis_sinusoidal_latlon(h::Integer, v::Integer,
                                 line::Integer, sample::Integer;
                                 pixels::Integer=MODIS_DEFAULT_PIXELS)
    _validate_modis_tile(h, v)
    _validate_modis_pixels(pixels)
    1 <= line <= pixels || error("MODIS line must be in 1:$pixels, got $line")
    1 <= sample <= pixels || error("MODIS sample must be in 1:$pixels, got $sample")

    radius = MODIS_SINUSOIDAL_RADIUS
    tile = modis_tile_width()
    pixel = tile / pixels

    x = -pi * radius + h * tile + (sample - 0.5) * pixel
    y = pi * radius / 2 - v * tile - (line - 0.5) * pixel

    lat_rad = y / radius
    lon_rad = x / (radius * cos(lat_rad))

    Float32(rad2deg(lat_rad)), Float32(rad2deg(lon_rad))
end

"""
    generate_modis_tile_geolocation(h, v; pixels=2400)

Generate `latitude` and `longitude` arrays for a MODLAND sinusoidal tile from
the official projection definition. Arrays are indexed as `[line, sample]`.
"""
function generate_modis_tile_geolocation(h::Integer, v::Integer;
                                         pixels::Integer=MODIS_DEFAULT_PIXELS)
    _validate_modis_tile(h, v)
    _validate_modis_pixels(pixels)

    lat = Matrix{Float32}(undef, pixels, pixels)
    lon = Matrix{Float32}(undef, pixels, pixels)

    radius = MODIS_SINUSOIDAL_RADIUS
    tile = modis_tile_width()
    pixel = tile / pixels
    x0 = -pi * radius + h * tile
    y0 = pi * radius / 2 - v * tile

    @inbounds for line in 1:pixels
        y = y0 - (line - 0.5) * pixel
        lat_rad = y / radius
        lat_deg = Float32(rad2deg(lat_rad))
        cos_lat = cos(lat_rad)
        for sample in 1:pixels
            x = x0 + (sample - 0.5) * pixel
            lat[line, sample] = lat_deg
            lon[line, sample] = Float32(rad2deg(x / (radius * cos_lat)))
        end
    end

    lat, lon
end

function modis_tile_cache_path(cache_dir::AbstractString, h::Integer, v::Integer;
                               pixels::Integer=MODIS_DEFAULT_PIXELS)
    _validate_modis_tile(h, v)
    _validate_modis_pixels(pixels)
    joinpath(String(cache_dir), @sprintf("h%02dv%02d_%dpx.nc", h, v, pixels))
end

"""
    write_modis_tile_geolocation(path, h, v; pixels=2400, overwrite=false)

Generate and write one MODIS sinusoidal tile geolocation file. The output is a
small per-tile NetCDF file with `latitude(y, x)` and `longitude(y, x)`.
"""
function write_modis_tile_geolocation(path::AbstractString, h::Integer, v::Integer;
                                      pixels::Integer=MODIS_DEFAULT_PIXELS,
                                      overwrite::Bool=false)
    _validate_modis_tile(h, v)
    _validate_modis_pixels(pixels)

    outpath = String(path)
    if isfile(outpath) && !overwrite
        return outpath
    end

    outdir = dirname(outpath)
    !isempty(outdir) && mkpath(outdir)
    lat, lon = generate_modis_tile_geolocation(h, v; pixels=pixels)

    tmp = outpath * ".tmp-$(getpid())"
    ds = Dataset(tmp, "c")
    try
        defDim(ds, "y", pixels)
        defDim(ds, "x", pixels)

        vlat = defVar(ds, "latitude", Float32, ("y", "x"),
                      deflatelevel=4, fillvalue=Float32(NaN))
        vlon = defVar(ds, "longitude", Float32, ("y", "x"),
                      deflatelevel=4, fillvalue=Float32(NaN))
        vlat[:, :] = lat
        vlon[:, :] = lon

        ds.attrib["title"] = "Generated MODIS sinusoidal tile geolocation"
        ds.attrib["modis_tile_h"] = h
        ds.attrib["modis_tile_v"] = v
        ds.attrib["pixels_per_tile"] = pixels
        ds.attrib["sphere_radius_m"] = MODIS_SINUSOIDAL_RADIUS
        ds.attrib["source"] = "Generated from MODLAND sinusoidal grid parameters"
        ds.attrib["source_url"] = "https://modis-land.gsfc.nasa.gov/GCTP.html"
    finally
        close(ds)
    end

    mv(tmp, outpath; force=overwrite)
    outpath
end

function read_modis_tile_geolocation(path::AbstractString)
    ds = Dataset(String(path))
    try
        ds["latitude"][:, :], ds["longitude"][:, :]
    finally
        close(ds)
    end
end

"""
    load_or_generate_modis_tile_geolocation(h, v; cache_dir=default_modis_cache_dir(), pixels=2400)

Read a per-tile MODIS geolocation cache file, generating it first if needed.
Returns `(latitude, longitude)` arrays indexed as `[line, sample]`.
"""
function load_or_generate_modis_tile_geolocation(h::Integer, v::Integer;
                                                 cache_dir::AbstractString=default_modis_cache_dir(),
                                                 pixels::Integer=MODIS_DEFAULT_PIXELS)
    path = modis_tile_cache_path(cache_dir, h, v; pixels=pixels)
    if !isfile(path)
        write_modis_tile_geolocation(path, h, v; pixels=pixels)
    end
    read_modis_tile_geolocation(path)
end

function center_geolocation_provider(config::DataSourceConfig;
                                     geo_table::Union{Nothing,String}=nothing,
                                     geo_cache::Union{Nothing,String}=nothing,
                                     geo_provider::Union{Symbol,String}=:auto)
    provider = Symbol(geo_provider)

    if provider == :auto
        if geo_table !== nothing
            return ModisMonolithicLUTGeolocation(geo_table)
        elseif haskey(config.basic, "lat") && haskey(config.basic, "lon")
            return BasicVariableGeolocation()
        else
            pixels = Int(get(config.options, "modis_pixels", MODIS_DEFAULT_PIXELS))
            cache_dir = geo_cache === nothing ? default_modis_cache_dir(; pixels=pixels) : geo_cache
            return ModisSinusoidalGeolocation(cache_dir, pixels)
        end
    elseif provider in (:variables, :basic)
        return BasicVariableGeolocation()
    elseif provider in (:lut, :geo_table, :geotable)
        geo_table === nothing && error("geo_provider='$provider' requires geo_table")
        return ModisMonolithicLUTGeolocation(geo_table)
    elseif provider in (:modis, :modis_sinusoidal)
        pixels = Int(get(config.options, "modis_pixels", MODIS_DEFAULT_PIXELS))
        cache_dir = geo_cache === nothing ? default_modis_cache_dir(; pixels=pixels) : geo_cache
        return ModisSinusoidalGeolocation(cache_dir, pixels)
    end

    error("Unknown center geolocation provider '$geo_provider'. Use auto, variables, lut, or modis.")
end

function open_geolocation!(provider::AbstractCenterGeolocation)
    provider
end

function open_geolocation!(provider::ModisMonolithicLUTGeolocation)
    provider.dataset = Dataset(provider.path)
    provider
end

function close_geolocation!(provider::AbstractCenterGeolocation)
    nothing
end

function close_geolocation!(provider::ModisMonolithicLUTGeolocation)
    provider.dataset !== nothing && close(provider.dataset)
    provider.dataset = nothing
    nothing
end

function center_coordinates(provider::BasicVariableGeolocation,
                            filepath::String, config::DataSourceConfig)
    haskey(config.basic, "lat") && haskey(config.basic, "lon") ||
        error("center geolocation provider 'variables' requires basic.lat and basic.lon")
    read_array_from_file(filepath, config.basic["lat"]),
    read_array_from_file(filepath, config.basic["lon"])
end

function center_coordinates(provider::ModisMonolithicLUTGeolocation,
                            filepath::String, config::DataSourceConfig)
    provider.dataset === nothing && error("MODIS geolocation LUT is not open")
    h, v = parse_modis_tile(filepath)
    provider.dataset["latitude"][h + 1, v + 1, :, :],
    provider.dataset["longitude"][h + 1, v + 1, :, :]
end

function center_coordinates(provider::ModisSinusoidalGeolocation,
                            filepath::String, config::DataSourceConfig)
    h, v = try
        parse_modis_tile(filepath)
    catch e
        error("center gridding requires basic.lat/basic.lon, geo_table, or a MODIS-style filename containing hXXvYY for generated sinusoidal geolocation: $e")
    end
    load_or_generate_modis_tile_geolocation(h, v;
                                            cache_dir=provider.cache_dir,
                                            pixels=provider.pixels)
end
