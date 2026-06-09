# Zarr v2 daily emit -----------------------------------------------------------
#
# Additive, per-day Zarr v2 store written alongside the NetCDF output, readable
# unchanged by Zarr.jl, xarray/zarr, and zarrita.js. One store grows over time;
# each day is one chunk per array.
#
# Zarr is a SOFT (weak) dependency: the array-creating/writing functions
# (`init_zarr_store`, `write_zarr_day!`) live in the `SatelliteGriddingZarrExt`
# package extension and only become available after `using Zarr`. The pure
# helpers below (layout transform, time index, consolidated metadata, meta.json)
# need no Zarr and always load. The public `init_zarr_store`/`write_zarr_day!`
# here are thin delegators to the extension (mirrors `cuda_backend`/`metal_backend`).
#
# ORDERING CONTRACT (the whole risk of this feature) -----------------------------
# The in-memory grid here is `[lon, lat, var]` with latitude ASCENDING (index 1 =
# −89.5, south) and longitude in [−180,180) (index 1 = −179.5). The store is
# logical `(time, lat, lon)` with latitude DESCENDING (row 0 = +89.5, north) and
# longitude ascending.
#
# KEY FACT: Zarr.jl writes Julia arrays with their dimensions REVERSED in the
# on-disk `.zarray` metadata (with `order:"C"`). A Julia ZArray of shape
# `(NX, NY, NT)` is therefore read by xarray/zarrita as `(NT, NY, NX)` =
# `(time, lat, lon)` — exactly the logical layout we want. So the data arrays are
# CREATED in reversed Julia order `(lon, lat, time)`, and the ONLY per-slice
# transform is flipping latitude south-up → north-up (`to_zarr_slice`, a reverse
# along the lat axis); Zarr.jl's dimension reversal supplies the rest. Nothing is
# ever "fixed" in the browser. The ramp-probe test (test/test_zarr_ordering.jl +
# test/zarr_probe/) enforces this across Julia, Python, and the browser.

# Time axis is pre-allocated to a generous horizon; Zarr only materializes chunks
# that are written, so empty future/past days cost nothing. EPOCH/HORIZON MUST be
# identical across every run that shares one store (they define the time axis).
const ZARR_EPOCH = Date(2021, 1, 1)
const ZARR_HORIZON = Date(2031, 1, 1)
const ZARR_NT = Dates.value(ZARR_HORIZON - ZARR_EPOCH) + 1

# Dashboard variable registry (label/units/vmin/vmax/cmap/mask) for meta.json.
# Unknown vars fall back to `_default_registry_entry`.
const ZARR_VAR_REGISTRY = Dict{String,Dict{String,Any}}(
    "xco2"     => Dict("label" => "XCO2", "units" => "ppm", "vmin" => 406, "vmax" => 424, "cmap" => "rdylbu_r", "mask" => "none"),
    "xco2_unc" => Dict("label" => "XCO2 uncertainty", "units" => "ppm", "vmin" => 0.2, "vmax" => 1.1, "cmap" => "cividis", "mask" => "none"),
    "sif757"   => Dict("label" => "SIF 757 nm", "units" => "W m-2 sr-1 um-1", "vmin" => 0, "vmax" => 2.2, "cmap" => "viridis", "mask" => "land"),
    "sif771"   => Dict("label" => "SIF 771 nm", "units" => "W m-2 sr-1 um-1", "vmin" => 0, "vmax" => 1.8, "cmap" => "viridis", "mask" => "land"),
    "sifdaily" => Dict("label" => "Daily-mean SIF 757 nm", "units" => "W m-2 sr-1 um-1", "vmin" => 0, "vmax" => 1.6, "cmap" => "magma", "mask" => "land"),
)

_default_registry_entry(v) = Dict{String,Any}(
    "label" => v, "units" => _zarr_units(v), "vmin" => 0, "vmax" => 1,
    "cmap" => "viridis", "mask" => "none")

# Best-effort units for the Zarr array `units` attribute (display units live in meta.json).
function _zarr_units(v::AbstractString)
    lv = lowercase(v)
    startswith(lv, "xco2") && return "ppm"
    occursin("sif", lv) && return "W m-2 sr-1 um-1"
    return ""
end

"""
    to_zarr_slice(slice_lonlat) -> Matrix

Prepare one in-memory `[lon, lat]` slice (size `(NX, NY)`, latitude ascending /
south-up) for writing as a Zarr `(lon, lat)` Julia slice with latitude DESCENDING
(north-up). This is the single ordering guard: it flips latitude (`reverse` along
dim 2). It does NOT transpose — the data array is created in reversed Julia order
`(lon, lat, time)` so Zarr.jl's dimension reversal makes xarray/zarrita read it as
`(time, lat, lon)`. Result is `(NX, NY)` with column 0 (lat) at +89.5.
"""
to_zarr_slice(slice_lonlat::AbstractMatrix) = reverse(slice_lonlat; dims = 2)

# Absolute, date-based time index (1-based) → identical across runs and re-runs,
# which is what makes per-day writes idempotent and date-parallel-safe.
zarr_time_index(date, epoch::Date = ZARR_EPOCH) = Dates.value(Date(date) - epoch) + 1

# True if a Zarr array `name` already exists under the store at `path`.
_zarr_array_exists(path, name) = isfile(joinpath(path, name, ".zarray"))

"""
    init_zarr_store(path, grid_spec, grid_vars; epoch, horizon, compute_std) -> ZGroup

Open-or-create the Zarr v2 store at `path` and ensure it has the arrays for every
variable in `grid_vars`. Additive: an existing store is opened and only the arrays
it lacks are created, so several gridding runs (e.g. SIF and XCO2) can write their
own variables into one combined store.

Requires `using Zarr` (provided by the `SatelliteGriddingZarrExt` extension).
"""
function init_zarr_store(args...; kwargs...)
    ext = Base.get_extension(@__MODULE__, :SatelliteGriddingZarrExt)
    ext === nothing &&
        error("Zarr daily emit requires Zarr.jl. Run `using Zarr` to load the SatelliteGriddingZarrExt extension (add Zarr to your environment first).")
    # invokelatest so it works even when Zarr was loaded at runtime (newer world age),
    # e.g. bin/grid.jl's `@eval using Zarr` when --zarrOut is set.
    Base.invokelatest(ext.init_zarr_store, args...; kwargs...)
end

"""
    write_zarr_day!(g, date, grid_data, grid_std, grid_weights, grid_vars, compute_std; epoch)

Write one day's gridded fields into the store `g` (the group returned by
[`init_zarr_store`](@ref)) at the absolute time index for `date`. `grid_data`,
`grid_std` are `[lon, lat, var]`; `grid_weights` is `[lon, lat]` (fractional
footprint weights). Writes exactly one chunk per array: the transposed/flipped
mean, the per-cell count (`Int16`, weights rounded; empty → −1), and — when
`compute_std` — the standard error of the mean `std / √n`. Empty cells stay
NaN / −1; zeros are never written.

Requires `using Zarr` (provided by the `SatelliteGriddingZarrExt` extension).
"""
function write_zarr_day!(args...; kwargs...)
    ext = Base.get_extension(@__MODULE__, :SatelliteGriddingZarrExt)
    ext === nothing &&
        error("Zarr daily emit requires Zarr.jl. Run `using Zarr` to load the SatelliteGriddingZarrExt extension (add Zarr to your environment first).")
    Base.invokelatest(ext.write_zarr_day!, args...; kwargs...)
end

"""
    consolidate_zarr!(path)

Regenerate consolidated metadata (`.zmetadata`) so the browser fetches all store
metadata in one request. Pure Julia (no Zarr needed): walk the store and embed
every `.zgroup`, `.zarray`, and `.zattrs`.
"""
function consolidate_zarr!(path::AbstractString)
    meta = Dict{String,Any}()
    for (root, _, files) in walkdir(path)
        for f in files
            if f in (".zgroup", ".zarray", ".zattrs")
                rel = replace(relpath(joinpath(root, f), path), '\\' => '/')
                meta[rel] = JSON.parsefile(joinpath(root, f))
            end
        end
    end
    open(joinpath(path, ".zmetadata"), "w") do io
        JSON.print(io, Dict("zarr_consolidated_format" => 1, "metadata" => meta))
    end
    return path
end

"""
    write_zarr_meta_json(webdir, grid_spec, grid_vars; epoch, registry)

Emit/merge `meta.json` (grid definition + per-variable dashboard registry) into
`webdir`. Additive: an existing `meta.json` is read and its `variables` are merged,
so the SIF and XCO2 runs both contribute their variables to one registry. Pure
Julia (no Zarr needed).
"""
function write_zarr_meta_json(webdir::AbstractString, grid_spec::RectangularGridSpec,
                              grid_vars; epoch::Date = ZARR_EPOCH,
                              registry = ZARR_VAR_REGISTRY)
    isdir(webdir) || mkpath(webdir)
    path = joinpath(webdir, "meta.json")

    variables = Dict{String,Any}()
    if isfile(path)
        prev = JSON.parsefile(path)
        haskey(prev, "variables") && merge!(variables, Dict{String,Any}(prev["variables"]))
    end
    for v in keys(grid_vars)
        variables[v] = get(registry, v, _default_registry_entry(v))
    end

    doc = Dict(
        "grid" => Dict("nx" => length(grid_spec.lon), "ny" => length(grid_spec.lat),
                       "lat0" => grid_spec.lat_max, "lon0" => grid_spec.lon_min,
                       "dLat" => grid_spec.dlat, "dLon" => grid_spec.dlon),
        "epoch" => string(epoch),
        "variables" => variables,
    )
    open(path, "w") do io
        JSON.print(io, doc)
    end
    return path
end
