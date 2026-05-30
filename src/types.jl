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

# The output-grid geometry types (`AbstractGridSpec`, `RectangularGridSpec`,
# `GridSpec` alias) and their dispatch interface live in `grid_spec.jl`. The
# top-level `SatelliteGridding.jl` includes that file immediately after this
# one so the configuration / metadata records above stay independent of the
# grid-extensibility surface.

"""
    TimeSpec

Specification of the temporal gridding parameters.

Two temporal modes are supported:

- **Tiled (default):** snapshots are spaced by `time_step`; each snapshot averages a
  forward window `[d, d + time_step·oversample_temporal − 1 day]`. With
  `oversample_temporal == 1` the windows tile with no overlap.
- **Rolling mean:** set `sample_step` and `window_halfwidth_days`. Snapshots are spaced
  by `sample_step` and each averages a *centered* window `[c − N, c + N]` (N =
  `window_halfwidth_days`). When `2N+1 > sample_step` the windows overlap, smoothing the
  temporal evolution. Each calendar day is gridded only once and reused across
  overlapping windows (in-memory day cache). `time_step`/`oversample_temporal` are
  ignored in this mode.

# Fields
- `start_date`, `stop_date`: Date range of snapshot centers
- `time_step`: Tiled-mode bin size (`Dates.Day` or `Dates.Month`)
- `oversample_temporal`: Tiled-mode window multiplier (>1 gives moving-average-like behavior)
- `sample_step`: Rolling-mode snapshot spacing (`Dates.Day`), or `nothing` for tiled mode
- `window_halfwidth_days`: Rolling-mode half-width N (days), or `nothing` for tiled mode
"""
struct TimeSpec
    start_date::DateTime
    stop_date::DateTime
    time_step::Union{Dates.Day,Dates.Month}
    oversample_temporal::Float32
    sample_step::Union{Nothing,Dates.Day}
    window_halfwidth_days::Union{Nothing,Int}
end

"""
    TimeSpec(start_date, stop_date, time_step; oversample_temporal=1.0f0,
             sample_step=nothing, window_halfwidth_days=nothing)

Construct a `TimeSpec`. Pass both `sample_step` and `window_halfwidth_days` to enable
rolling-mean mode; pass neither for the default tiled mode.
"""
function TimeSpec(start_date::DateTime, stop_date::DateTime,
                  time_step::Union{Dates.Day,Dates.Month};
                  oversample_temporal::Float32=1.0f0,
                  sample_step::Union{Nothing,Dates.Day}=nothing,
                  window_halfwidth_days::Union{Nothing,Int}=nothing)
    (sample_step === nothing) == (window_halfwidth_days === nothing) ||
        error("Rolling mean requires both sample_step and window_halfwidth_days, or neither.")
    if window_halfwidth_days !== nothing && window_halfwidth_days < 0
        error("window_halfwidth_days must be ≥ 0; got $window_halfwidth_days.")
    end
    TimeSpec(start_date, stop_date, time_step, oversample_temporal,
             sample_step, window_halfwidth_days)
end

"""
    is_rolling(ts::TimeSpec) -> Bool

True when `ts` is configured for centered rolling-mean gridding.
"""
is_rolling(ts::TimeSpec) = ts.sample_step !== nothing

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
