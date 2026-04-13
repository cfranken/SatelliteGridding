# SatelliteGridding.jl

A Julia package for gridding satellite Level-2 data onto regular rectangular grids
with footprint-aware oversampling.

## Overview

Satellite instruments observe through irregularly shaped footprints — typically
quadrilaterals defined by 4 corner coordinates. These footprints can overlap multiple
output grid cells. `SatelliteGridding.jl` handles this by subdividing each footprint
into n×n sub-pixels and using "point-in-box" binning to distribute observations across
the grid.

### Key Features

- **Footprint oversampling**: Subdivide footprints into n×n sub-pixels for proper
  spatial overlap handling
- **TOML configuration**: Declarative data source definitions with filter expressions
- **Date-based file discovery**: Automatic file matching using YYYY/MM/DD/DOY placeholders
- **Quality filtering**: Flexible filter syntax (`"< 80"`, `"> 0.5"`, `"== 1"`, `"1600 < x < 2200"`)
- **Multiple backends**: Sequential (Welford), KA CPU, or KA GPU (CUDA) execution
- **Standard deviation**: Optional per-cell standard deviation via Welford's algorithm
- **Center-coordinate gridding**: MODIS-style gridding where each pixel maps to one cell
- **Vegetation indices**: Built-in EVI, NDVI, NIRv, NDWI computation for MODIS data

### Supported Instruments

- **TROPOMI**: SIF, NO₂, CH₄, CO (S5P Level-2)
- **OCO-2/3**: SIF, XCO₂ (Level-2)
- **MODIS**: Surface reflectance (center-coordinate gridding)
- Any Level-2 product with corner coordinates in NetCDF format

## Prerequisites

**Julia ≥ 1.9** is required. If you don't have Julia installed:

1. Download from [julialang.org/downloads](https://julialang.org/downloads/)
2. Follow the platform-specific install instructions
3. Verify with `julia --version`

## Installation

From the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/cfranken/SatelliteGridding.git")
```

This installs the package and all dependencies automatically.

For development (editable local clone):

```julia
Pkg.develop(url="https://github.com/cfranken/SatelliteGridding.git")
```

## Quick Start — Julia API

```julia
using SatelliteGridding, Dates

# 1. Load a TOML config that describes the satellite data source
config = load_config("examples/tropomi_sif.toml")

# 2. Define the output grid (global, 0.5° resolution)
grid_spec = GridSpec(
    lat_min = -90f0,  lat_max = 90f0,
    lon_min = -180f0, lon_max = 180f0,
    dlat = 0.5f0, dlon = 0.5f0
)

# 3. Define the time range (8-day composites)
time_spec = TimeSpec(
    DateTime("2019-01-01"),
    DateTime("2019-12-31"),
    Dates.Day(8)
)

# 4. Run the gridding
grid_l2(config, grid_spec, time_spec;
        outfile = "tropomi_sif_2019.nc")
```

### Reading the output

```julia
using NCDatasets

ds = Dataset("tropomi_sif_2019.nc")
sif     = ds["sif_743"][:, :, 1]   # gridded SIF, first time step
weights = ds["n"][:, :, 1]          # observation count per cell
close(ds)
```

### Key options for [`grid_l2`](@ref)

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `n_oversample` | `nothing` (auto) | Sub-pixel factor. `nothing` = auto-compute from footprint/grid ratio |
| `compute_std` | `false` | Compute per-cell standard deviation (sequential backend only) |
| `outfile` | `"gridded_output.nc"` | Output NetCDF file path |
| `backend` | `nothing` (sequential) | `CPU()` for KA parallel, `CUDABackend()` for GPU |

## Quick Start — Command Line

The CLI has two commands: **`l2`** (footprint oversampling) and **`center`** (MODIS-style).

```bash
# TROPOMI SIF — 8-day composites at 0.5°
julia --project=. bin/grid.jl l2 \
    --config SatelliteGridding/examples/tropomi_sif.toml \
    --dLat 0.5 --dLon 0.5 --dDays 8 \
    --startDate 2019-01-01 --stopDate 2019-12-31 \
    -o tropomi_sif_2019.nc

# OCO-2 XCO₂ — monthly composites at 2°
julia --project=. bin/grid.jl l2 \
    --config SatelliteGridding/examples/oco2_xco2.toml \
    --dLat 2.0 --dLon 2.0 --monthly --dDays 1 \
    --startDate 2015-01-01 --stopDate 2020-12-31 \
    -o oco2_xco2_monthly.nc

# With parallel KA CPU backend
julia --project=. bin/grid.jl l2 \
    --config SatelliteGridding/examples/tropomi_no2.toml \
    --backend cpu --dLat 0.25 --dLon 0.25 \
    -o tropomi_no2.nc
```

See the [CLI Reference](@ref) for all options and the [GPU Acceleration](@ref) page
for CUDA setup.

## Backends

| Backend | CLI flag | Julia API | Description |
|:--------|:---------|:----------|:------------|
| Sequential | `--backend sequential` | default (`backend=nothing`) | Welford's online mean. Supports `--compSTD`. |
| KA CPU | `--backend cpu` | `backend=CPU()` | Parallel sort + subpixels. Sum-based accumulation. |
| KA CUDA | `--backend cuda` | `backend=CUDABackend()` | Full GPU pipeline. Requires `using CUDA`. |

## Next Steps

- [Tutorial](@ref "Tutorial: Gridding TROPOMI SIF Data"): Step-by-step walkthrough with TROPOMI SIF data
- [Configuration](@ref "Configuration Reference"): TOML config reference and filter expressions
- [Algorithm](@ref): How footprint oversampling works under the hood
- [CLI Reference](@ref): All command-line flags
- [GPU Acceleration](@ref): CUDA setup and performance guide
- [API Reference](@ref): Full function documentation
