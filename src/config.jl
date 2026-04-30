"""
    load_config(config_path::String) -> DataSourceConfig

Load a configuration file (TOML or JSON) and return a `DataSourceConfig`.

The file must contain a `[grid]` section plus `filePattern` and `folder`.
Footprint-aware L2 gridding usually uses `[basic]` entries for `lat_bnd` and
`lon_bnd`. Circular footprints can alternatively use center `lat`/`lon` plus a
`radius` variable or scalar `[circle] radius`. Center-coordinate gridding can use
`[basic]` `lat` and `lon`, a geolocation lookup table, or generated MODIS
sinusoidal geolocation. Filters are defined in an optional `[filter]` section
using intuitive string expressions:

# Example TOML
```toml
filePattern = "S5P_PAL__L2B_SIF____YYYYMMDD*.nc"
folder = "/path/to/data/"

[basic]
lat = "PRODUCT/latitude"
lon = "PRODUCT/longitude"
lat_bnd = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds"
lon_bnd = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds"

[grid]
sif_743 = "PRODUCT/SIF_743"
sif_735 = "PRODUCT/SIF_735"

[filter]
"PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle" = "< 80"
"PRODUCT/qa_value" = "> 0.5"
"PRODUCT/methane" = "1600 < x < 2200"
"PRODUCT/mode" = "== 1"
```

The legacy JSON format with separate `filter_gt`/`filter_lt`/`filter_eq` sections
is still supported for backward compatibility.
"""
function load_config(config_path::String)::DataSourceConfig
    if endswith(config_path, ".json")
        d = JSON.parsefile(config_path)
    else
        d = TOML.parsefile(config_path)
    end

    # Required sections. Center-coordinate gridding can derive geolocation from
    # a separate lookup table, so `[basic]` is optional there.
    haskey(d, "grid") || error("Config missing required 'grid' section: $config_path")

    basic = Dict{String,String}()
    if haskey(d, "basic")
        basic = Dict{String,String}(k => string(v) for (k, v) in d["basic"])
    end

    # Preserve key order for grid variables (band order matters for MODIS)
    grid_vars = OrderedDict{String,String}()
    for (k, v) in d["grid"]
        grid_vars[k] = string(v)
    end

    # Parse filters
    filters = FilterRule[]

    # New unified [filter] section with string expressions
    if haskey(d, "filter")
        for (var, expr) in d["filter"]
            push!(filters, _parse_filter_expr(var, expr))
        end
    end

    # Legacy format: separate filter_gt / filter_lt / filter_eq sections
    for (key, value) in get(d, "filter_gt", Dict())
        push!(filters, FilterRule(string(key), :gt, Float64(value)))
    end
    for (key, value) in get(d, "filter_lt", Dict())
        push!(filters, FilterRule(string(key), :lt, Float64(value)))
    end
    for (key, value) in get(d, "filter_eq", Dict())
        push!(filters, FilterRule(string(key), :eq, Float64(value)))
    end

    file_pattern = string(get(d, "filePattern", ""))
    folder = string(get(d, "folder", ""))
    !isempty(file_pattern) || error("Config missing required 'filePattern': $config_path")
    !isempty(folder) || error("Config missing required 'folder': $config_path")

    options = Dict{String,Any}()
    for section in ("options", "center", "modis", "circle")
        if haskey(d, section)
            for (k, v) in d[section]
                options[string(k)] = v
            end
        end
    end
    _warn_unknown_options(options)

    DataSourceConfig(basic, grid_vars, filters, file_pattern, folder, options)
end

function _warn_unknown_options(options::Dict{String,Any})
    known = Set([
        "scale_factor", "add_offset", "fill_value", "valid_min", "valid_max",
        "transpose_data", "min_count", "min_nir_reflectance", "modis_pixels",
        "vegetation_red", "vegetation_nir", "vegetation_blue", "vegetation_swir",
        "radius", "radius_unit",
    ])
    for key in keys(options)
        if !(key in known)
            @warn "Unknown config option; it will be ignored unless custom code reads it" option=key
        end
    end
    nothing
end

"""
    _parse_filter_expr(variable, expr) -> FilterRule

Parse a filter expression string into a `FilterRule`.

Supported formats:
- `"< 80"` → less than
- `"> 0.5"` → greater than
- `"== 1"` → equality
- `"1600 < x < 2200"` → range (between, exclusive)
"""
function _parse_filter_expr(variable::String, expr)::FilterRule
    s = strip(string(expr))

    # Range: "1600 < x < 2200" or "1600 <x< 2200"
    m = match(r"^(-?[\d.eE+\-]+)\s*<\s*\w+\s*<\s*(-?[\d.eE+\-]+)$", s)
    if m !== nothing
        return FilterRule(variable, :between, parse(Float64, m[1]), parse(Float64, m[2]))
    end

    # Less than: "< 80"
    m = match(r"^<\s*(-?[\d.eE+\-]+)$", s)
    if m !== nothing
        return FilterRule(variable, :lt, parse(Float64, m[1]))
    end

    # Greater than: "> 0.5"
    m = match(r"^>\s*(-?[\d.eE+\-]+)$", s)
    if m !== nothing
        return FilterRule(variable, :gt, parse(Float64, m[1]))
    end

    # Equality: "== 1"
    m = match(r"^==\s*(-?[\d.eE+\-]+)$", s)
    if m !== nothing
        return FilterRule(variable, :eq, parse(Float64, m[1]))
    end

    # Plain number (from TOML integer/float value, not a string) — treat as equality
    m = match(r"^(-?[\d.eE+\-]+)$", s)
    if m !== nothing
        return FilterRule(variable, :eq, parse(Float64, m[1]))
    end

    error("Cannot parse filter expression for '$variable': \"$expr\"\n" *
          "Expected: \"< 80\", \"> 0.5\", \"== 1\", or \"1600 < x < 2200\"")
end
