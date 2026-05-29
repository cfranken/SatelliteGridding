#!/usr/bin/env julia
#
# Grid TROPOMI L2 SIF onto a GEOS/GCHP cubed-sphere grid (library usage).
#
# The data-source config (variable paths, filters, file pattern) is grid-agnostic
# and shared with the regular lat/lon examples — only the `grid_spec` differs.
# Here we use a C24 grid for a quick demo; set Nc=360 for a GCHP-style C360 run.
#
#   julia --project=. examples/tropomi_sif_cubedsphere.jl
#
# Output is a GCHP-style NetCDF (dims nf/Ydim/Xdim, 2D lons/lats + corner arrays)
# that opens directly in Panoply, xarray, and GCHP tooling.

using Dates
using SatelliteGridding

config = load_config(joinpath(@__DIR__, "tropomi_sif.toml"))

# GMAOCubedSphereDefinition() is the default: GEOS/GCHP equal-distance gnomonic,
# GEOS-native panel ordering, -10° longitude shift. Use Nc=360 for C360.
grid_spec = CubedSphereGridSpec(; Nc = 24, definition = GMAOCubedSphereDefinition())

time_spec = TimeSpec(
    DateTime("2020-07-01"),
    DateTime("2020-07-16"),
    Dates.Day(16),
)

# Cubed-sphere gridding currently runs on the sequential CPU path only.
grid(config, grid_spec, time_spec, SubpixelGridding();
     keep_going = true,
     outfile = "tropomi_sif_c24_2020-07-01.nc")
