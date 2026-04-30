"""
    grid(config, grid_spec, time_spec; method=SubpixelGridding(), kwargs...)
    grid(config, grid_spec, time_spec, method::AbstractGriddingMethod; kwargs...)

Dispatch-based library entry point for gridding. The CLI still exposes explicit
`l2` and `center` commands, while library callers can choose a gridding method
with a first-class type such as `SubpixelGridding`,
`CircularFootprintGridding`, or `CenterPointGridding`.
"""
function grid(config::DataSourceConfig, grid_spec::GridSpec,
              time_spec::TimeSpec; method::AbstractGriddingMethod=SubpixelGridding(),
              kwargs...)
    grid(config, grid_spec, time_spec, method; kwargs...)
end

function grid(config::DataSourceConfig, grid_spec::GridSpec,
              time_spec::TimeSpec, method::SubpixelGridding; kwargs...)
    grid_l2(config, grid_spec, time_spec;
            footprint_method=method,
            kwargs...)
end

function grid(config::DataSourceConfig, grid_spec::GridSpec,
              time_spec::TimeSpec, method::CircularFootprintGridding; kwargs...)
    grid_l2(config, grid_spec, time_spec;
            footprint_method=method,
            kwargs...)
end

function grid(config::DataSourceConfig, grid_spec::GridSpec,
              time_spec::TimeSpec, ::CenterPointGridding; kwargs...)
    grid_center(config, grid_spec, time_spec; kwargs...)
end

function grid(config::DataSourceConfig, grid_spec::GridSpec,
              time_spec::TimeSpec, ::ExactIntersectionGridding; kwargs...)
    error("ExactIntersectionGridding is not implemented yet. Use SubpixelGridding for approximate footprint overlap or CenterPointGridding for fast center binning.")
end
