#!/usr/bin/env julia

using Dates
using SatelliteGridding

config = load_config(joinpath(@__DIR__, "gosat_sif_center_radius.toml"))

grid_spec = GridSpec(
    lat_min = -90.0f0,
    lat_max = 90.0f0,
    lon_min = -180.0f0,
    lon_max = 180.0f0,
    dlat = 1.0f0,
    dlon = 1.0f0,
)

time_spec = TimeSpec(
    DateTime("2010-01-01"),
    DateTime("2010-12-31"),
    Dates.Day(8),
)

backend_name = get(ENV, "SATGRID_BACKEND", "sequential")
backend = resolve_backend(backend_name)

grid(config, grid_spec, time_spec,
     CircularFootprintGridding(n_oversample=80);
     backend=backend,
     outfile="gosat_sif_2010_center_radius.nc")
