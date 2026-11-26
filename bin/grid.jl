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

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."); io="--quiet" in ARGS ? devnull : stderr)
end

using Dates
using SatelliteGridding

function main(argv=ARGS; default_progress=SatelliteGridding._default_progress())
    if isempty(argv) || argv[1] in ["-h", "--help"]
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
    command_line = join(["julia", PROGRAM_FILE, argv...], " ")
    command = first(argv)
    args = if command == "l2"
        parse_l2_args(argv[2:end])
    elseif command == "center"
        parse_center_args(argv[2:end])
    else
        error("Unknown command: $command. Use 'l2' or 'center'.")
    end
    config = load_config(args["config"])
    grid_spec = args_to_grid_spec(args)
    time_spec = args_to_time_spec(args)
    run = (outfile=args["outFile"], keep_going=args["keepGoing"],
           progress=!args["noProgress"] && default_progress, verbose=!args["quiet"],
           command=command_line)

    if command == "l2"
        geometry_precision = SatelliteGridding._cli_precision(args["geometry-precision"])
        accumulator_precision = SatelliteGridding._cli_precision(args["accumulator-precision"])
        method = SatelliteGridding._cli_footprint(args)
        backend = resolve_backend(args["backend"]; grid_spec, geometry_precision, accumulator_precision)

        zarr_out = isempty(args["zarrOut"]) ? nothing : args["zarrOut"]
        zarr_meta = isempty(args["zarrMeta"]) ? nothing : args["zarrMeta"]
        # Loading Zarr triggers the SatelliteGriddingZarrExt extension that backs
        # init_zarr_store / write_zarr_day!. Only load it when a Zarr store is requested.
        if zarr_out !== nothing
            SatelliteGridding._validate_zarr_output(grid_spec,
                SatelliteGridding._output_options(config.options))
            @eval using Zarr
        end

        grid(config, grid_spec, time_spec, method;
             compute_std=args["compSTD"], zarr_out, zarr_meta,
             zarr_epoch=Date(args["zarrEpoch"]), backend, geometry_precision, accumulator_precision,
             run...)
    else
        geo = isempty(args["geoTable"]) ? nothing : args["geoTable"]
        geo_cache = isempty(args["geoCache"]) ? nothing : args["geoCache"]

        grid(config, grid_spec, time_spec, CenterPointGridding();
             geo_table=geo, geo_cache, geo_provider=args["geoProvider"],
             veg_indices=args["vegIndices"], run...)
    end

    if !args["quiet"]
        elapsed = time() - t_start
        println("Total time: $(round(elapsed, digits=1)) seconds")
    end
    nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
