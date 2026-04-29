"""
    SatelliteGridding

A Julia package for gridding satellite Level-2 data onto regular rectangular grids
with footprint-aware oversampling.

Satellite instruments have irregularly shaped footprints (quadrilaterals defined by
4 corner coordinates) that can overlap multiple output grid cells. Rather than computing
exact geometric intersections, this package subdivides each footprint into n×n sub-pixels
and uses "point-in-box" binning. The mean per grid cell is computed incrementally using
Welford's online algorithm for numerical stability.

# Main functions
- [`grid_l2`](@ref): Grid Level-2 data with footprint oversampling
- [`grid_center`](@ref): Grid data using center coordinates only (e.g., MODIS)
- [`load_config`](@ref): Load a TOML/JSON data source configuration

# Example
```julia
using SatelliteGridding

config = load_config("tropomi_sif.toml")
grid_spec = GridSpec(lat_min=-90f0, lat_max=90f0, lon_min=-180f0, lon_max=180f0,
                     dlat=1f0, dlon=1f0)
time_spec = TimeSpec(DateTime("2019-01-01"), DateTime("2019-12-31"),
                     Dates.Day(8), 1.0f0)
grid_l2(config, grid_spec, time_spec; outfile="output.nc")
```
"""
module SatelliteGridding

using Dates
using Printf
using Statistics
using TOML
using JSON
using OrderedCollections
using Glob
using NCDatasets
using ArgParse
using ProgressMeter
using KernelAbstractions
using Atomix: @atomic
using Libdl

include("types.jl")
include("config.jl")
include("filepatterns.jl")
include("system_netcdf.jl")
include("ncio.jl")
include("geolocation.jl")
include("filters.jl")
include("oversampling.jl")
include("geometry.jl")
include("averaging.jl")
include("vegetation.jl")
include("gridder.jl")
include("methods.jl")
include("cli.jl")

export GridSpec, DataSourceConfig, TimeSpec, FilterRule
export AbstractGriddingMethod, SubpixelGridding, CenterPointGridding,
       ExactIntersectionGridding
export load_config, find_files
export read_nc_variable, read_nc_attribute, read_variable_from_file, read_array_from_file,
       create_output_dataset
export read_system_netcdf_variable
export AbstractCenterGeolocation, BasicVariableGeolocation,
       ModisMonolithicLUTGeolocation, ModisSinusoidalGeolocation
export default_cache_dir, default_modis_cache_dir, parse_modis_tile,
       modis_sinusoidal_latlon, generate_modis_tile_geolocation,
       modis_tile_cache_path, write_modis_tile_geolocation,
       read_modis_tile_geolocation, load_or_generate_modis_tile_geolocation
export apply_filters
export compute_subpixels!, floor_indices!, compute_n_oversample
export polygon_area, exact_footprint_weights, footprint_weights
export sort_corners_ccw, sort_corners_ccw!, sort_corners_ccw_ka!
export compute_footprint_indices_ka!, scatter_accumulate!, scatter_accumulate_ka!
export accumulate_footprint!, accumulate_batch!, accumulate_center!
export finalize_mean!
export grid, grid_l2, grid_center
export compute_evi, compute_ndvi, compute_nirv, compute_ndwi
export parse_l2_args, parse_center_args, args_to_grid_spec, args_to_time_spec

end # module
