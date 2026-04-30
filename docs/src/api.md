# API Reference

## Module

```@docs
SatelliteGridding
```

## Main Functions

```@docs
grid
grid_l2
grid_center
load_config
```

## Backends

```@docs
supported_backend_names
backend_help_lines
resolve_backend
cuda_backend
metal_backend
```

## Types

```@docs
GridSpec
TimeSpec
DataSourceConfig
FilterRule
AbstractGriddingMethod
SubpixelGridding
CircularFootprintGridding
CenterPointGridding
ExactIntersectionGridding
```

## Oversampling

```@docs
sort_corners_ccw
sort_corners_ccw!
sort_corners_ccw_ka!
compute_subpixels!
floor_indices!
compute_n_oversample
compute_footprint_indices_ka!
compute_circular_footprint_indices_ka!
```

## Geometry

```@docs
polygon_area
exact_footprint_weights
circle_footprint_weights
footprint_weights
```

## Accumulation

```@docs
accumulate_footprint!
accumulate_circular_footprint!
accumulate_batch!
accumulate_circular_batch!
accumulate_center!
scatter_accumulate!
scatter_accumulate_ka!
scatter_accumulate_circular!
scatter_accumulate_circular_ka!
finalize_mean!
```

## I/O

```@docs
find_files
read_nc_variable
read_nc_attribute
read_variable_from_file
read_array_from_file
read_system_netcdf_variable
create_output_dataset
```

## Geolocation

```@docs
default_cache_dir
default_modis_cache_dir
parse_modis_tile
modis_sinusoidal_latlon
generate_modis_tile_geolocation
write_modis_tile_geolocation
load_or_generate_modis_tile_geolocation
```

## Filters

```@docs
apply_filters
```

## Vegetation Indices

```@docs
compute_evi
compute_ndvi
compute_nirv
compute_ndwi
```

## CLI

```@docs
parse_l2_args
parse_center_args
args_to_grid_spec
args_to_time_spec
```
