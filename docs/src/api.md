# API Reference

## Main Functions

```@docs
grid_l2
grid_center
load_config
```

## Types

```@docs
GridSpec
TimeSpec
DataSourceConfig
FilterRule
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
```

## Accumulation

```@docs
accumulate_footprint!
accumulate_batch!
accumulate_center!
scatter_accumulate!
scatter_accumulate_ka!
finalize_mean!
```

## I/O

```@docs
find_files
read_nc_variable
read_nc_attribute
create_output_dataset
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
