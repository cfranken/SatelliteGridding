#!/usr/bin/env julia

# CLI entry point for SatelliteGridding
#
# Usage:
#   julia --project=. bin/grid.jl l2 --config examples/tropomi_sif.toml --dLat 0.5 --dLon 0.5 -o output.nc
#   julia --project=. bin/grid.jl l2 --config examples/tropomi_sif.toml --backend cpu -o output.nc
#   julia --project=. bin/grid.jl center --config examples/modis_reflectance.toml --geoProvider modis -o modis_out.nc
#
# For help:
#   julia --project=. bin/grid.jl l2 --help
#   julia --project=. bin/grid.jl center --help

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using SatelliteGridding

function main()
    if length(ARGS) < 1 || ARGS[1] in ["-h", "--help"]
        println("SatelliteGridding CLI v0.1.0")
        println()
        println("Usage:")
        println("  julia --project=. bin/grid.jl <command> [options]")
        println()
        println("Commands:")
        println("  l2      Grid Level-2 data with footprint oversampling")
        println("  center  Grid data using center coordinates (MODIS-style)")
        println()
        println("Backends:")
        for line in backend_help_lines()
            println("  --backend ", line)
        end
        println()
        println("Run 'julia --project=. bin/grid.jl <command> --help' for command-specific options.")
        return
    end

    t_start = time()
    command = popfirst!(ARGS)

    if command == "l2"
        args = parse_l2_args(ARGS)
        config = load_config(args["config"])
        grid_spec = args_to_grid_spec(args)
        time_spec = args_to_time_spec(args)
        n_over = args["nOversample"] > 0 ? args["nOversample"] : nothing
        backend = resolve_backend(args["backend"])
        footprint = lowercase(args["footprint"])
        method = if footprint in ("quad", "quadrilateral", "subpixel")
            SubpixelGridding(n_over)
        elseif footprint in ("circle", "circular")
            CircularFootprintGridding(n_over)
        else
            error("Unknown footprint geometry: $(args["footprint"]). Use quad or circle.")
        end

        grid(config, grid_spec, time_spec, method;
             compute_std=args["compSTD"],
             outfile=args["outFile"],
             backend=backend,
             keep_going=args["keepGoing"])

    elseif command == "center"
        args = parse_center_args(ARGS)
        config = load_config(args["config"])
        grid_spec = args_to_grid_spec(args)
        time_spec = args_to_time_spec(args)
        geo = isempty(args["geoTable"]) ? nothing : args["geoTable"]
        geo_cache = isempty(args["geoCache"]) ? nothing : args["geoCache"]

        grid(config, grid_spec, time_spec, CenterPointGridding();
             geo_table=geo,
             geo_cache=geo_cache,
             geo_provider=args["geoProvider"],
             veg_indices=args["vegIndices"],
             outfile=args["outFile"],
             keep_going=args["keepGoing"])
    else
        error("Unknown command: $command. Use 'l2' or 'center'.")
    end

    elapsed = time() - t_start
    println("Total time: $(round(elapsed, digits=1)) seconds")
end

main()
