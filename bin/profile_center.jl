using SatelliteGridding
using SatelliteGridding: _center_values, _center_grid_indices,
                          apply_center_filters, center_coordinates,
                          accumulate_center!, sort_files_for_provider
using Profile, Glob

usage() = error("usage: julia --project=. bin/profile_center.jl <config.toml> <DOY-dir> [n_files=10]")

length(ARGS) >= 2 || usage()
config_path = ARGS[1]
doy_dir     = ARGS[2]
n_files     = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10

config     = load_config(config_path)
grid_spec  = GridSpec(lat_min=-90.0f0, lat_max=90.0f0,
                     lon_min=-180.0f0, lon_max=180.0f0,
                     dlat=1.0f0, dlon=1.0f0)
provider   = ModisSinusoidalGeolocation(
    expanduser("~/.cache/SatelliteGridding/modis/sinusoidal_2400px_v1"), 2400)
varpaths   = collect(values(config.grid_vars))

files_all  = sort_files_for_provider(provider, glob("MCD43A4.*.hdf", doy_dir))
files      = files_all[1:min(n_files, length(files_all))]

println("config:    ", config_path)
println("doy_dir:   ", doy_dir)
println("vars:      ", length(varpaths))
println("files:     ", length(files), " (of ", length(files_all), ")")
println()

# Warm up Julia compile + first-touch overhead with one ignored file
let f = files[1]
    SatelliteGridding.read_arrays_from_file(f, varpaths)
    lat, lon = center_coordinates(provider, f, config)
    apply_center_filters(lat, lon, grid_spec)
end

totals = Dict("read"=>0.0, "geo"=>0.0, "filter"=>0.0, "idx"=>0.0,
              "inner"=>0.0, "veg"=>0.0, "wall"=>0.0)

println(rpad("file", 60), rpad("read", 9), rpad("geo", 8),
        rpad("filter", 9), rpad("idx", 8), rpad("inner", 9),
        rpad("wall", 9), "valid")

for f in files
    GC.gc()
    t_wall = time()

    t_read = @elapsed arrays = SatelliteGridding.read_arrays_from_file(f, varpaths)
    t_geo  = @elapsed (lat_in, lon_in) = center_coordinates(provider, f, config)
    t_flt  = @elapsed idx = apply_center_filters(lat_in, lon_in, grid_spec)
    t_idx  = @elapsed (lat_idx, lon_idx) = _center_grid_indices(idx, lat_in, lon_in, grid_spec)
    # _center_values reads again; bypass it for the inner-loop measurement by
    # mimicking its work without re-reading.
    n_vars = length(varpaths)
    n_p    = length(idx)
    values = zeros(Float32, n_p, n_vars + 4)
    valid  = trues(n_p)
    scale  = Float32(get(config.options, "scale_factor", 1.0))
    offset = Float32(get(config.options, "add_offset", 0.0))
    fill_v = get(config.options, "fill_value", nothing)
    vmin   = get(config.options, "valid_min", nothing)
    vmax   = get(config.options, "valid_max", nothing)
    transp = Bool(get(config.options, "transpose_data", false))

    t_inner = @elapsed begin
        for co in 1:n_vars
            SatelliteGridding._fill_center_column!(values, valid, idx,
                arrays[co], co, scale, offset, fill_v, vmin, vmax, transp)
        end
    end

    t_total_wall = time() - t_wall

    totals["read"]   += t_read
    totals["geo"]    += t_geo
    totals["filter"] += t_flt
    totals["idx"]    += t_idx
    totals["inner"]  += t_inner
    totals["wall"]   += t_total_wall

    println(rpad(basename(f), 60),
            rpad(string(round(Int, t_read*1000),  "ms"), 9),
            rpad(string(round(Int, t_geo*1000),   "ms"), 8),
            rpad(string(round(Int, t_flt*1000),   "ms"), 9),
            rpad(string(round(Int, t_idx*1000),   "ms"), 8),
            rpad(string(round(Int, t_inner*1000), "ms"), 9),
            rpad(string(round(Int, t_total_wall*1000), "ms"), 9),
            count(valid))
end

println()
n = length(files)
println("== averages over $n files ==")
for k in ("read","geo","filter","idx","inner","wall")
    println("  ", rpad(k, 8), round(Int, totals[k]/n*1000), " ms")
end
println("  ", rpad("sum", 8),
        round(Int, (totals["read"]+totals["geo"]+totals["filter"]+
                    totals["idx"]+totals["inner"])/n*1000), " ms (parts)")
println("  files/sec: ", round(n / totals["wall"], digits=2))

println()
println("== profile of the inner per-pixel loop (last file's scale) ==")
Profile.clear()
let f = files[end]
    arrays = SatelliteGridding.read_arrays_from_file(f, varpaths)
    lat_in, lon_in = center_coordinates(provider, f, config)
    idx = apply_center_filters(lat_in, lon_in, grid_spec)
    n_p = length(idx)
    values = zeros(Float32, n_p, length(varpaths) + 4)
    valid  = trues(n_p)
    scale  = Float32(get(config.options, "scale_factor", 1.0))
    offset = Float32(get(config.options, "add_offset", 0.0))
    fill_v = get(config.options, "fill_value", nothing)
    vmin   = get(config.options, "valid_min", nothing)
    vmax   = get(config.options, "valid_max", nothing)
    transp = Bool(get(config.options, "transpose_data", false))

    @profile for _ in 1:5
        for co in 1:length(varpaths)
            data = arrays[co]
            @inbounds for (i, ci) in enumerate(idx)
                raw = SatelliteGridding._center_data_value(data, ci, transp)
                if SatelliteGridding._invalid_center_value(raw, fill_v, vmin, vmax)
                    valid[i] = false
                    values[i, co] = Float32(NaN)
                else
                    values[i, co] = Float32(raw) * scale + offset
                end
            end
        end
    end
end
Profile.print(format=:flat, sortedby=:count, mincount=20, maxdepth=12)
