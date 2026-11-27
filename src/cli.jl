function _add_common_args!(s; spacing, start_date, stop_date, step_days, config_help, outfile_help)
    @add_arg_table! s begin
        "--config", "-c"
            help = config_help
            arg_type = String
            required = true
        "--outFile", "-o"
            help = outfile_help
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
            default = spacing
        "--dLon"
            help = "Longitude resolution (degrees)"
            arg_type = Float32
            default = spacing
        "--startDate"
            help = "Start date (YYYY-MM-DD)"
            arg_type = String
            default = start_date
        "--stopDate"
            help = "Stop date (YYYY-MM-DD)"
            arg_type = String
            default = stop_date
        "--dDays"
            help = "Time step in days (or months if --monthly is set)"
            arg_type = Int64
            default = step_days
        "--keepGoing"
            help = "Continue after per-file processing errors and report failures at the end"
            action = :store_true
        "--noProgress"
            help = "Don't draw progress bars. By default they're drawn only when stderr is a terminal"
            action = :store_true
        "--quiet"
            help = "Suppress progress bars and run status messages"
            action = :store_true
    end
    s
end

function _parse_grid_args(args, settings)
    parsed = parse_args(args, settings)
    parsed["noProgress"] |= parsed["quiet"]
    parsed
end

function _parse_coarsen_args(args)
    s = ArgParseSettings(description="Coarsen a NESTED HEALPix NetCDF archive")
    @add_arg_table! s begin
        "--input"
            help = "Input HEALPix NetCDF file"
            required = true
        "--output"
            help = "Output HEALPix NetCDF file"
            required = true
        "--levels"
            help = "Number of parent levels to coarsen (default: 1)"
            arg_type = Int
        "--nside"
            help = "Target nside, a power of two smaller than the input nside"
            arg_type = Int
        "--legacy-std"
            help = "Untagged <name>_std pairs: error (default), by_name (population standard deviations), or none (independent means). Tagged archives use their explicit links."
            default = "error"
            range_tester = x -> x in ("error", "by_name", "none")
    end
    parse_args(args, s)
end

"""
    parse_l2_args(args=ARGS) -> Dict

Parse command-line arguments for L2 gridding. Matches the interface of the
original `gridL2_Dates.jl` script.
"""
function parse_l2_args(args=ARGS)
    s = ArgParseSettings(description="Grid Level-2 satellite data with footprint oversampling")
    _add_common_args!(s; spacing=1.0f0, start_date="2018-03-07", stop_date="2018-10-31", step_days=8,
                      config_help="TOML or JSON configuration file defining the data source",
                      outfile_help="Output NetCDF filename; empty string disables NetCDF (Zarr-only, needs --zarrOut)")

    @add_arg_table! s begin
        "--compSTD"
            help = "Compute standard deviation within each grid cell"
            action = :store_true
        "--gridType"
            help = "Output grid geometry: rect (regular lat/lon), cs (GEOS/GCHP cubed sphere), or healpix (NESTED HEALPix)"
            arg_type = String
            default = "rect"
        "--Nc"
            help = "Cubed-sphere cells per panel edge (e.g. 360 for C360); required when --gridType cs"
            arg_type = Int64
            default = 0
        "--nside"
            help = "HEALPix resolution, a positive power of two required when --gridType healpix"
            arg_type = Int64
            default = 0
        "--csConvention"
            help = "Cubed-sphere convention: gmao (GEOS/GCHP, -10° offset) or equiangular"
            arg_type = String
            default = "gmao"
        "--oversample_temporal"
            help = "Temporal oversampling factor (averaging window = oversample_temporal × dDays)"
            arg_type = Float32
            default = 1.0f0
        "--sample_step"
            help = "Rolling mean: snapshot spacing in days. >0 enables centered rolling-mean mode and overrides dDays/oversample_temporal for scheduling."
            arg_type = Int64
            default = 0
        "--window_days"
            help = "Rolling mean: averaging half-width N in days; each snapshot averages the centered window [center-N, center+N]. Used only when --sample_step > 0."
            arg_type = Int64
            default = 0
        "--nOversample"
            help = "Sub-pixel factor for footprint oversampling (0 or omitted: automatic; 1–$N_OVERRIDE_MAX: fixed)"
            arg_type = Int64
        "--footprint"
            help = "Footprint geometry: quad for corner quadrilaterals, circle for circular footprints"
            arg_type = String
            default = "quad"
        "--backend"
            help = "Compute backend: sequential, cpu, cuda, or metal"
            arg_type = String
            default = "sequential"
        "--geometry-precision"
            help = "Spherical KA geometry precision: float32 or float64 (default: float32 on Metal, float64 elsewhere)"
            arg_type = String
            range_tester = x -> lowercase(x) in ("float32", "float64")
        "--accumulator-precision"
            help = "Accumulator precision: float32 or float64. Sets grid precision on the sequential backend. Sequential rolling requires Float64 daily accumulation and rejects explicit float32."
            arg_type = String
            range_tester = x -> lowercase(x) in ("float32", "float64")
        "--zarrOut"
            help = "Also write a per-day Zarr v2 store at this path (additive; one chunk per variable per day, alongside the NetCDF). Requires Zarr.jl. Rectangular grids only."
            arg_type = String
            default = ""
        "--zarrMeta"
            help = "Directory to write the dashboard meta.json (grid definition + variable registry). Used only with --zarrOut."
            arg_type = String
            default = ""
        "--zarrEpoch"
            help = "Zarr time-axis epoch (YYYY-MM-DD). Must be identical across all runs that share one store."
            arg_type = String
            default = "2021-01-01"
    end
    parsed = _parse_grid_args(args, s)
    if parsed["sample_step"] > 0 && _cli_sequential_backend(parsed) &&
       _cli_precision(parsed["accumulator-precision"]) === Float32
        throw(ArgumentError("--accumulator-precision float32 is unsupported for sequential rolling, which accumulates daily moments in Float64"))
    end
    n = parsed["nOversample"]
    parsed["nOversample"] = _validate_n_oversample(n == 0 ? nothing : n,
        "--nOversample must be 0 (automatic) or an integer in 1:$N_OVERRIDE_MAX")
    parsed
end

"""
    parse_center_args(args=ARGS) -> Dict

Parse command-line arguments for center-coordinate gridding (MODIS-style).
"""
function parse_center_args(args=ARGS)
    s = ArgParseSettings(description="Grid satellite data using center coordinates (no footprint bounds)")
    _add_common_args!(s; spacing=0.5f0, start_date="2018-01-01", stop_date="2018-12-31", step_days=1,
                      config_help="TOML or JSON configuration file", outfile_help="Output NetCDF filename")

    @add_arg_table! s begin
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
    end
    _parse_grid_args(args, s)
end

"""
    args_to_grid_spec(args::Dict) -> AbstractGridSpec

Convert parsed CLI arguments to a grid spec. Returns a `RectangularGridSpec` for
`--gridType rect` (default), a `CubedSphereGridSpec` for `--gridType cs`
(GEOS/GCHP convention by default, requires `--Nc`), or a `HEALPixGridSpec`
for `--gridType healpix` with `--nside`.
"""
function args_to_grid_spec(args::Dict)::AbstractGridSpec
    gtype = lowercase(strip(get(args, "gridType", "rect")))
    precision = _cli_sequential_backend(args) ? _cli_precision(get(args, "accumulator-precision", nothing)) : nothing
    T = something(precision, Float32)
    if gtype in ("rect", "rectangular", "latlon", "lonlat")
        convert_bound = something(precision, identity)
        return GridSpec(lat_min=convert_bound(args["latMin"]), lat_max=convert_bound(args["latMax"]),
                        lon_min=convert_bound(args["lonMin"]), lon_max=convert_bound(args["lonMax"]),
                        dlat=convert_bound(args["dLat"]), dlon=convert_bound(args["dLon"]))
    elseif gtype in ("cs", "cube", "cubed", "cubedsphere", "cubed_sphere")
        Nc = get(args, "Nc", 0)
        Nc > 0 || error("--gridType cs requires --Nc > 0 (e.g. --Nc 360 for C360)")
        conv = lowercase(strip(get(args, "csConvention", "gmao")))
        definition = if conv in ("gmao", "geos", "gchp")
            GMAOCubedSphereDefinition()
        elseif conv in ("equiangular", "equiangular_gnomonic", "synthetic")
            EquiangularCubedSphereDefinition()
        else
            error("Unknown --csConvention '$conv'. Use gmao or equiangular.")
        end
        return CubedSphereGridSpec(; Nc=Nc, definition, T)
    elseif gtype == "healpix"
        return HEALPixGridSpec(; nside=get(args, "nside", 0), T)
    else
        error("Unknown --gridType '$gtype'. Use rect, cs, or healpix.")
    end
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
    sample_step = get(args, "sample_step", 0)
    if sample_step > 0
        return TimeSpec(start_date, stop_date, time_step;
                        oversample_temporal=oversample,
                        sample_step=Dates.Day(sample_step),
                        window_halfwidth_days=get(args, "window_days", 0))
    end
    TimeSpec(start_date, stop_date, time_step; oversample_temporal=oversample)
end


_cli_precision(::Nothing) = nothing
_cli_precision(arg::AbstractString) = lowercase(arg) == "float32" ? Float32 : Float64
_cli_sequential_backend(args) = _is_sequential_backend(_normalize_backend_name(get(args, "backend", "sequential")))

function _cli_footprint(args)
    footprint = lowercase(args["footprint"])
    n = args["nOversample"]
    footprint in ("quad", "quadrilateral", "subpixel") && return SubpixelGridding(n)
    footprint in ("circle", "circular") && return CircularFootprintGridding(n)
    error("Unknown footprint geometry: $(args["footprint"]). Use quad or circle.")
end
