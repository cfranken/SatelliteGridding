#!/usr/bin/env julia
#
# Ramp-probe writer — the mandatory ordering gate (spec §4.1).
#
# Builds an in-memory [lon, lat] south-up grid engineered so that, AFTER the real
# writer's transpose+flip (`to_zarr_latlon` inside `write_zarr_day!`), the logical
# store cell (iy, ix) — 0-based, iy=0 NORTH, ix=0 lon=−180 — holds iy*1000 + ix.
# A single read then pins the entire index→coordinate mapping; any transpose,
# latitude flip, longitude roll, or C/F-order mistake fails loudly in check.py /
# check.mjs.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using SatelliteGridding
using Zarr
using Dates
using OrderedCollections

path = get(ARGS, 1, "probe.zarr")
isdir(path) && rm(path; recursive = true)

gs = GridSpec()                       # 1° global: 180 lat × 360 lon
gv = OrderedDict("probe" => "probe")

# G_mem[ix1, iy1] (1-based, [lon,lat], lat ascending/south-up). After
# reverse(permutedims(·,(2,1)); dims=1) this lands at store (iy=iy1... see derivation):
#   store a[iy, ix] (0-based) = iy*1000 + ix.
G_mem = Float32[(180 - iy1) * 1000 + (ix1 - 1) for ix1 in 1:360, iy1 in 1:180]
grid_data = reshape(G_mem, 360, 180, 1)
grid_w = ones(Float32, 360, 180)               # all populated → nothing masked
grid_std = zeros(Float32, 360, 180, 1)

g = init_zarr_store(path, gs, gv; compute_std = false)
write_zarr_day!(g, Date(2021, 1, 1), grid_data, grid_std, grid_w, gv, false)
consolidate_zarr!(path)
println("probe written to $path  (day index 1)")
