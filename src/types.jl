"""
    FilterRule

A single filter condition for a NetCDF variable.

# Fields
- `variable`: Path to the NetCDF variable (may include group separators `/`)
- `op`: Operation — `:lt`, `:gt`, `:eq`, `:between`
- `lo`: Lower bound (threshold for `:gt`/`:lt`/`:eq`, low end for `:between`)
- `hi`: Upper bound (only used for `:between`)

# Config syntax (TOML `[filter]` section)
```toml
[filter]
"solar_zenith_angle" = "< 80"         # less than
"qa_value" = "> 0.5"                  # greater than
"mode" = "== 1"                       # equality
"methane" = "1600 < x < 2200"         # range (between)
```
"""
struct FilterRule
    variable::String
    op::Symbol
    lo::Float64
    hi::Float64
end

FilterRule(variable::String, op::Symbol, value::Real) =
    FilterRule(variable, op, Float64(value), NaN)

"""
    DataSourceConfig

Configuration loaded from a TOML/JSON file that defines a satellite data source.

# Fields
- `basic`: Maps internal keys (`"lat"`, `"lon"`, `"lat_bnd"`, `"lon_bnd"`, `"radius"`) to variable paths in the NetCDF files
- `grid_vars`: Ordered mapping of output variable names to input variable paths (all will be gridded)
- `filters`: Quality filter rules parsed from `[filter]` section
- `file_pattern`: Glob pattern with YYYY/MM/DD/DOY placeholders for finding input files
- `folder`: Root folder for input data (may also contain YYYY/MM/DD placeholders)
- `options`: Optional processing settings from `[options]`, `[center]`, or `[modis]`
"""
struct DataSourceConfig
    basic::Dict{String,String}
    grid_vars::OrderedDict{String,String}
    filters::Vector{FilterRule}
    file_pattern::String
    folder::String
    options::Dict{String,Any}
end

DataSourceConfig(basic::Dict{String,String}, grid_vars::OrderedDict{String,String},
                 filters::Vector{FilterRule}, file_pattern::String, folder::String) =
    DataSourceConfig(basic, grid_vars, filters, file_pattern, folder, Dict{String,Any}())

"""
    GridSpec{T<:AbstractFloat}

Specification of the output grid geometry.

# Fields
- `lat_min`, `lat_max`: Latitude bounds (degrees)
- `lon_min`, `lon_max`: Longitude bounds (degrees)
- `dlat`, `dlon`: Grid cell size (degrees)
- `lat`, `lon`: Vectors of cell center coordinates
"""
struct GridSpec{T<:AbstractFloat}
    lat_min::T
    lat_max::T
    lon_min::T
    lon_max::T
    dlat::T
    dlon::T
    lat::Vector{T}
    lon::Vector{T}
end

"""
    GridSpec(; lat_min=-90f0, lat_max=90f0, lon_min=-180f0, lon_max=180f0, dlat=1f0, dlon=1f0)

Construct a `GridSpec` with cell centers computed from bounds and resolution.
"""
function GridSpec(; lat_min::T=-90.0f0, lat_max::T=90.0f0,
                   lon_min::T=-180.0f0, lon_max::T=180.0f0,
                   dlat::T=1.0f0, dlon::T=1.0f0) where {T<:AbstractFloat}
    eps = dlat / 100
    lat = collect(lat_min + dlat / 2:dlat:lat_max - dlat / 2 + eps)
    lon = collect(lon_min + dlon / 2:dlon:lon_max - dlon / 2 + eps)
    GridSpec{T}(lat_min, lat_max, lon_min, lon_max, dlat, dlon, lat, lon)
end

"""
    TimeSpec

Specification of the temporal gridding parameters.

# Fields
- `start_date`, `stop_date`: Date range for processing
- `time_step`: Temporal bin size (`Dates.Day` or `Dates.Month`)
- `oversample_temporal`: Multiplier for the actual averaging window (>1 gives moving-average-like behavior)
"""
struct TimeSpec
    start_date::DateTime
    stop_date::DateTime
    time_step::Union{Dates.Day,Dates.Month}
    oversample_temporal::Float32
end

"""
    TimeSpec(start_date, stop_date, time_step; oversample_temporal=1.0f0)

Construct a `TimeSpec`.
"""
function TimeSpec(start_date::DateTime, stop_date::DateTime,
                  time_step::Union{Dates.Day,Dates.Month};
                  oversample_temporal::Float32=1.0f0)
    TimeSpec(start_date, stop_date, time_step, oversample_temporal)
end

"""
    AbstractGriddingMethod

Dispatch marker for choosing how observations are mapped onto the output grid.
Use [`SubpixelGridding`](@ref) for footprint-aware oversampling,
[`CircularFootprintGridding`](@ref) for circular or near-circular footprints,
[`CenterPointGridding`](@ref) for fast center-coordinate binning, and
[`ExactIntersectionGridding`](@ref) as the future exact geometry hook.
"""
abstract type AbstractGriddingMethod end

"""
    SubpixelGridding(; n_oversample=nothing)

Footprint-aware gridding by sampling an `n × n` grid of subpoints inside each
quadrilateral footprint. `nothing` keeps the existing automatic oversampling
heuristic.
"""
struct SubpixelGridding <: AbstractGriddingMethod
    n_oversample::Union{Nothing,Int}
end

SubpixelGridding(; n_oversample::Union{Nothing,Int}=nothing) =
    SubpixelGridding(n_oversample)

"""
    CircularFootprintGridding(; n_oversample=nothing)

Footprint-aware gridding for circular or near-circular footprints described by
a center coordinate and four bounding coordinates. The implementation samples an
`n × n` grid over the footprint bounding box and keeps points inside the inferred
circle or ellipse. `nothing` keeps the existing automatic oversampling heuristic.
"""
struct CircularFootprintGridding <: AbstractGriddingMethod
    n_oversample::Union{Nothing,Int}
end

CircularFootprintGridding(; n_oversample::Union{Nothing,Int}=nothing) =
    CircularFootprintGridding(n_oversample)

"""
    CenterPointGridding()

Fast center-coordinate gridding. Each observation or pixel contributes to the
single output cell containing its center coordinate.
"""
struct CenterPointGridding <: AbstractGriddingMethod end

"""
    ExactIntersectionGridding()

Reserved method for exact footprint/grid-cell intersection gridding. This is not
implemented yet, but the type keeps the public API open for an exact geometry
backend without changing callers.
"""
struct ExactIntersectionGridding <: AbstractGriddingMethod end
