"""
    parse_l2_args(args=ARGS) -> Dict

Parse command-line arguments for L2 gridding. Matches the interface of the
original `gridL2_Dates.jl` script.
"""
function parse_l2_args(args=ARGS)
    s = ArgParseSettings(description="Grid Level-2 satellite data onto regular grids with footprint oversampling")

    @add_arg_table! s begin
        "--config", "-c"
            help = "TOML or JSON configuration file defining the data source"
            arg_type = String
            required = true
        "--outFile", "-o"
            help = "Output NetCDF filename"
            arg_type = String
            default = "gridded_output.nc"
        "--monthly"
            help = "Use time-steps in months instead of days"
            action = :store_true
        "--compSTD"
            help = "Compute standard deviation within each grid cell"
            action = :store_true
        "--latMin"
            help = "Lower latitude bound"
            arg_type = Float32
            default = -90.0f0
        "--latMax"
            help = "Upper latitude bound"
            arg_type = Float32
            default = 90.0f0
        "--lonMin"
            help = "Lower longitude bound"
            arg_type = Float32
            default = -180.0f0
        "--lonMax"
            help = "Upper longitude bound"
            arg_type = Float32
            default = 180.0f0
        "--dLat"
            help = "Latitude resolution (degrees)"
            arg_type = Float32
            default = 1.0f0
        "--dLon"
            help = "Longitude resolution (degrees)"
            arg_type = Float32
            default = 1.0f0
        "--startDate"
            help = "Start date (YYYY-MM-DD)"
            arg_type = String
            default = "2018-03-07"
        "--stopDate"
            help = "Stop date (YYYY-MM-DD)"
            arg_type = String
            default = "2018-10-31"
        "--dDays"
            help = "Time step in days (or months if --monthly is set)"
            arg_type = Int64
            default = 8
        "--oversample_temporal"
            help = "Temporal oversampling factor (averaging window = oversample_temporal × dDays)"
            arg_type = Float32
            default = 1.0f0
        "--nOversample"
            help = "Sub-pixel factor for footprint oversampling (default: auto-compute)"
            arg_type = Int64
            default = 0
        "--footprint"
            help = "Footprint geometry: quad for corner quadrilaterals, circle for circular footprints"
            arg_type = String
            default = "quad"
        "--backend"
            help = "Compute backend: sequential, cpu, cuda, or metal"
            arg_type = String
            default = "sequential"
        "--keepGoing"
            help = "Continue after per-file processing errors and report failures at the end"
            action = :store_true
    end
    parse_args(args, s)
end

"""
    parse_center_args(args=ARGS) -> Dict

Parse command-line arguments for center-coordinate gridding (MODIS-style).
"""
function parse_center_args(args=ARGS)
    s = ArgParseSettings(description="Grid satellite data using center coordinates (no footprint bounds)")

    @add_arg_table! s begin
        "--config", "-c"
            help = "TOML or JSON configuration file"
            arg_type = String
            required = true
        "--outFile", "-o"
            help = "Output NetCDF filename"
            arg_type = String
            default = "gridded_output.nc"
        "--monthly"
            help = "Use time-steps in months instead of days"
            action = :store_true
        "--latMin"
            help = "Lower latitude bound"
            arg_type = Float32
            default = -90.0f0
        "--latMax"
            help = "Upper latitude bound"
            arg_type = Float32
            default = 90.0f0
        "--lonMin"
            help = "Lower longitude bound"
            arg_type = Float32
            default = -180.0f0
        "--lonMax"
            help = "Upper longitude bound"
            arg_type = Float32
            default = 180.0f0
        "--dLat"
            help = "Latitude resolution (degrees)"
            arg_type = Float32
            default = 0.5f0
        "--dLon"
            help = "Longitude resolution (degrees)"
            arg_type = Float32
            default = 0.5f0
        "--startDate"
            help = "Start date (YYYY-MM-DD)"
            arg_type = String
            default = "2018-01-01"
        "--stopDate"
            help = "Stop date (YYYY-MM-DD)"
            arg_type = String
            default = "2018-12-31"
        "--dDays"
            help = "Time step in days (or months if --monthly is set)"
            arg_type = Int64
            default = 1
        "--geoTable"
            help = "Path to legacy monolithic geolocation lookup table (NetCDF)"
            arg_type = String
            default = ""
        "--geoCache"
            help = "Directory for generated per-tile MODIS sinusoidal geolocation cache"
            arg_type = String
            default = ""
        "--geoProvider"
            help = "Center geolocation provider: auto, variables, lut, or modis"
            arg_type = String
            default = "auto"
        "--vegIndices"
            help = "Compute vegetation indices (EVI, NDVI, NIRv, NDWI)"
            action = :store_true
        "--keepGoing"
            help = "Continue after per-file processing errors and report failures at the end"
            action = :store_true
    end
    parse_args(args, s)
end

"""
    args_to_grid_spec(args::Dict) -> GridSpec{Float32}

Convert parsed CLI arguments to a `GridSpec`.
"""
function args_to_grid_spec(args::Dict)::GridSpec{Float32}
    GridSpec(lat_min=args["latMin"], lat_max=args["latMax"],
             lon_min=args["lonMin"], lon_max=args["lonMax"],
             dlat=args["dLat"], dlon=args["dLon"])
end

"""
    args_to_time_spec(args::Dict) -> TimeSpec

Convert parsed CLI arguments to a `TimeSpec`.
"""
function args_to_time_spec(args::Dict)::TimeSpec
    start_date = DateTime(args["startDate"])
    stop_date = DateTime(args["stopDate"])
    time_step = get(args, "monthly", false) ? Dates.Month(args["dDays"]) : Dates.Day(args["dDays"])
    oversample = get(args, "oversample_temporal", 1.0f0)
    TimeSpec(start_date, stop_date, time_step; oversample_temporal=oversample)
end
