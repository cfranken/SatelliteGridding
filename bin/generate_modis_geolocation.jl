#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ArgParse
using SatelliteGridding

function parse_args_generate(args=ARGS)
    s = ArgParseSettings(description="Generate per-tile MODIS sinusoidal geolocation cache files")

    @add_arg_table! s begin
        "--tiles"
            help = "Comma-separated tile IDs, for example h08v04,h09v04"
            arg_type = String
            default = ""
        "--h"
            help = "Single horizontal MODIS tile index"
            arg_type = Int
            default = -1
        "--v"
            help = "Single vertical MODIS tile index"
            arg_type = Int
            default = -1
        "--range"
            help = "Generate a rectangular tile range using --hMin/--hMax/--vMin/--vMax"
            action = :store_true
        "--hMin"
            help = "Minimum horizontal tile for --range"
            arg_type = Int
            default = 0
        "--hMax"
            help = "Maximum horizontal tile for --range"
            arg_type = Int
            default = 35
        "--vMin"
            help = "Minimum vertical tile for --range"
            arg_type = Int
            default = 0
        "--vMax"
            help = "Maximum vertical tile for --range"
            arg_type = Int
            default = 17
        "--all"
            help = "Generate all 36 x 18 MODIS sinusoidal tiles"
            action = :store_true
        "--pixels"
            help = "Pixels per tile edge; 2400 for 500 m MCD43A4"
            arg_type = Int
            default = 2400
        "--cacheDir"
            help = "Output cache directory"
            arg_type = String
            default = ""
        "--force"
            help = "Overwrite existing tile files"
            action = :store_true
    end

    parse_args(args, s)
end

function requested_tiles(args)
    tiles = Tuple{Int,Int}[]

    if !isempty(args["tiles"])
        for tile in split(args["tiles"], ",")
            push!(tiles, parse_modis_tile(strip(tile)))
        end
    end

    if args["h"] >= 0 || args["v"] >= 0
        args["h"] >= 0 && args["v"] >= 0 ||
            error("Use both --h and --v for a single tile")
        push!(tiles, (args["h"], args["v"]))
    end

    if args["range"]
        for h in args["hMin"]:args["hMax"], v in args["vMin"]:args["vMax"]
            push!(tiles, (h, v))
        end
    end

    if args["all"]
        for h in 0:35, v in 0:17
            push!(tiles, (h, v))
        end
    end

    unique(tiles)
end

function main()
    args = parse_args_generate()
    tiles = requested_tiles(args)
    isempty(tiles) &&
        error("No tiles requested. Use --tiles h08v04, --h/--v, --range, or --all.")

    pixels = args["pixels"]
    cache_dir = isempty(args["cacheDir"]) ?
        default_modis_cache_dir(; pixels=pixels) : args["cacheDir"]

    println("Writing $(length(tiles)) MODIS geolocation tile(s) to $cache_dir")
    for (i, (h, v)) in enumerate(tiles)
        path = modis_tile_cache_path(cache_dir, h, v; pixels=pixels)
        println("[$i/$(length(tiles))] h$(lpad(h, 2, '0'))v$(lpad(v, 2, '0')) -> $path")
        write_modis_tile_geolocation(path, h, v;
                                     pixels=pixels,
                                     overwrite=args["force"])
    end
end

main()
