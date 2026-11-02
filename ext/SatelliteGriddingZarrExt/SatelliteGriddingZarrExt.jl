module SatelliteGriddingZarrExt

using SatelliteGridding
using Zarr
using Dates

# Pure helpers + constants that live in the parent (SatelliteGridding/src/zarrio.jl).
# `to_zarr_slice`, `zarr_time_index`, `RectangularGridSpec` are exported; the rest
# are internal and qualified via import.
import SatelliteGridding: RectangularGridSpec, to_zarr_slice, zarr_time_index,
                          ZARR_EPOCH, ZARR_HORIZON,
                          _zarr_units, _zarr_array_exists

# When Zarr is loaded, the parent's `init_zarr_store`/`write_zarr_day!` delegators
# call these implementations via `Base.get_extension(...).init_zarr_store(...)`.
# (Mirrors how SatelliteGriddingCUDAExt provides `cuda_backend`.)
#
# IMPORTANT (see src/zarrio.jl ORDERING CONTRACT): Zarr.jl reverses array
# dimensions in the on-disk `.zarray` metadata. Data arrays are therefore created
# in reversed Julia order `(NX, NY, NT)` = (lon, lat, time) so xarray/zarrita read
# them as `(NT, NY, NX)` = (time, lat, lon). Each day is written as `Z[:, :, ti]`.

"""
    init_zarr_store(path, grid_spec, grid_vars; epoch, horizon, compute_std) -> ZGroup

See `SatelliteGridding.init_zarr_store`. Open-or-create, additive: an existing store
is opened and only missing arrays are created. Data arrays are Julia `(NX, NY, NT)`
chunked `(NX, NY, 1)` (→ xarray/zarrita `(NT, NY, NX)` chunk `(1, NY, NX)`). For each
variable `v`: `v` (Float32, NaN fill), `v_n` (Int16, −1 fill), and — when
`compute_std` — `v_sem` (Float32, NaN fill), each with a CF `_ARRAY_DIMENSIONS`.
"""
function init_zarr_store(path::AbstractString, grid_spec::RectangularGridSpec,
                         grid_vars; epoch::Date = ZARR_EPOCH,
                         horizon::Date = ZARR_HORIZON, compute_std::Bool = true)
    NY = length(grid_spec.lat)
    NX = length(grid_spec.lon)
    NT = Dates.value(horizon - epoch) + 1
    comp = Zarr.BloscCompressor(cname = "zstd", clevel = 5, shuffle = 1)

    fresh = !(isdir(path) && isfile(joinpath(path, ".zgroup")))
    g = fresh ?
        zgroup(path; attrs = Dict("Conventions" => "CF-1.10",
                                  "title" => "OCO-2 daily gridded SIF & XCO2",
                                  "created_with" => "SatelliteGridding.jl")) :
        zopen(path, "w")

    if fresh
        # 1-D coordinates are unaffected by Zarr.jl's dimension reversal.
        t = zcreate(Float64, g, "time", NT; chunks = (NT,),
            attrs = Dict("_ARRAY_DIMENSIONS" => ["time"],
                         "units" => "days since $(epoch)",
                         "calendar" => "proleptic_gregorian"))
        t[:] = collect(0.0:(NT - 1))

        lat = zcreate(Float64, g, "lat", NY; chunks = (NY,),
            attrs = Dict("_ARRAY_DIMENSIONS" => ["lat"], "units" => "degrees_north",
                         "standard_name" => "latitude"))
        # DESCENDING, north-up: grid_spec.lat is ascending south-up, so reverse it
        # (exact — avoids float-range off-by-one at non-1° resolutions).
        lat[:] = Float64.(reverse(grid_spec.lat))

        lon = zcreate(Float64, g, "lon", NX; chunks = (NX,),
            attrs = Dict("_ARRAY_DIMENSIONS" => ["lon"], "units" => "degrees_east",
                         "standard_name" => "longitude"))
        lon[:] = collect(Float64.(grid_spec.lon))   # ASCENDING east
    end

    dims = ["time", "lat", "lon"]                    # logical (xarray/zarrita) order
    for v in keys(grid_vars)
        units = _zarr_units(v)
        if !_zarr_array_exists(path, v)
            zcreate(Float32, g, v, NX, NY, NT; chunks = (NX, NY, 1),
                compressor = comp, fill_value = NaN32,
                attrs = Dict("_ARRAY_DIMENSIONS" => dims, "units" => units))
        end
        if !_zarr_array_exists(path, v * "_n")
            zcreate(Int16, g, v * "_n", NX, NY, NT; chunks = (NX, NY, 1),
                compressor = comp, fill_value = Int16(-1),
                attrs = Dict("_ARRAY_DIMENSIONS" => dims, "units" => "count",
                             "long_name" => "Number of pixels in average"))
        end
        if compute_std && !_zarr_array_exists(path, v * "_sem")
            zcreate(Float32, g, v * "_sem", NX, NY, NT; chunks = (NX, NY, 1),
                compressor = comp, fill_value = NaN32,
                attrs = Dict("_ARRAY_DIMENSIONS" => dims, "units" => units,
                             "long_name" => "Standard error of the mean"))
        end
    end
    return g
end

"""
    write_zarr_day!(g, date, grid_data, grid_std, grid_weights, grid_vars, compute_std; epoch) -> ti

See `SatelliteGridding.write_zarr_day!`. Writes one `(NX, NY)` lat-flipped slice per
array at the absolute Julia time index `ti` (`Z[:, :, ti]`), i.e. exactly one chunk.
"""
function write_zarr_day!(g, date, grid_data::AbstractArray{<:Any,3},
                         grid_std::AbstractArray{<:Any,3}, grid_weights::AbstractMatrix,
                         grid_vars, compute_std::Bool; epoch::Date = ZARR_EPOCH)
    ti = zarr_time_index(date, epoch)
    # Bound against the store's ACTUAL time axis (sized from its own epoch/horizon),
    # not the default constant — `epoch` here must match the store's epoch.
    NT = size(g["time"])[1]
    1 <= ti <= NT ||
        error("date $date is outside the Zarr store time axis (time index $ti, NT=$NT); check the epoch matches the store, or widen its horizon.")

    wll = Float32.(to_zarr_slice(grid_weights))    # (NX, NY): lon × lat(north-up), fractional weights
    NX, NY = size(wll)
    # Validate against the store's actual spatial dims (works for any rectangular
    # resolution, not just 1°); also catches a transpose / wrong grid_spec.
    sNX, sNY, _ = size(g[first(keys(grid_vars))])
    (NX, NY) == (sNX, sNY) ||
        error("lat-flipped slice $((NX, NY)) does not match the store grid $((sNX, sNY)) — transpose bug or grid_spec mismatch?")
    empty = wll .< 1f-10

    # Int16 per-cell count: rounded weight where populated, −1 fill where empty.
    cnt = Array{Int16}(undef, NX, NY)
    @inbounds for i in eachindex(cnt)
        w = wll[i]
        cnt[i] = w < 1f-10 ? Int16(-1) : round(Int16, clamp(w, 0f0, 32767f0))
    end

    for (z, v) in enumerate(keys(grid_vars))
        val = Float32.(to_zarr_slice(@view grid_data[:, :, z]))
        @inbounds for i in eachindex(val)
            (empty[i] || !isfinite(val[i])) && (val[i] = NaN32)
        end
        g[v][:, :, ti] = val
        g[v * "_n"][:, :, ti] = cnt

        if compute_std
            sem = Float32.(to_zarr_slice(@view grid_std[:, :, z]))
            @inbounds for i in eachindex(sem)
                if empty[i]
                    sem[i] = NaN32
                else
                    sem[i] = sem[i] / sqrt(max(wll[i], eps(Float32)))
                    isfinite(sem[i]) || (sem[i] = NaN32)
                end
            end
            g[v * "_sem"][:, :, ti] = sem
        end
    end
    return ti
end

end # module
