# SatelliteGridding.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://cfranken.github.io/SatelliteGridding/dev/)

A Julia package for gridding satellite Level-2 data onto regular lat/lon grids with footprint-aware oversampling.

Satellite instruments observe through irregularly shaped footprints (quadrilaterals defined by 4 corner coordinates) that can overlap multiple grid cells. This package subdivides each footprint into n×n sub-pixels and distributes the observation across all grid cells the sub-pixels fall in. It supports any Level-2 satellite product with corner coordinates in NetCDF format — TROPOMI SIF/NO₂/CH₄/CO, OCO-2/3 SIF/XCO₂, GOSAT, and more.

## Prerequisites

**Julia ≥ 1.9** is required.

If you don't have Julia installed:

1. Download from [julialang.org/downloads](https://julialang.org/downloads/)
2. Follow the platform-specific install instructions
3. Verify: `julia --version`

MODIS HDF4/HDF-EOS2 files may require a system `libnetcdf` built with HDF4
support. The package first tries normal `NCDatasets.jl` reads and then falls
back to the system NetCDF C library for flat MODIS variables.

## Installation

From the Julia REPL (`julia` or press `]` for Pkg mode):

```julia
using Pkg
Pkg.add(url="https://github.com/cfranken/SatelliteGridding.git")
```

This installs the package and all its dependencies (NCDatasets, KernelAbstractions, ArgParse, etc.) automatically.

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

# 3. Define the time range (8-day composites for a full year)
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

### Key options for `grid_l2`

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `n_oversample` | `nothing` (auto) | Sub-pixel factor per footprint side. `nothing` = auto-compute from footprint/grid ratio |
| `compute_std` | `false` | Also compute per-cell standard deviation (sequential backend only) |
| `outfile` | `"gridded_output.nc"` | Output NetCDF file path |
| `backend` | `nothing` (sequential) | Compute backend — `nothing`, `CPU()`, or `CUDABackend()` (see [Backends](#backends)) |

### Center-coordinate gridding (MODIS)

For instruments without footprint bounds, use `grid_center`. MODIS sinusoidal
tile geolocation can be generated and cached per tile, so a large global lat/lon
lookup table is not required:

```julia
config = load_config("examples/modis_reflectance.toml")
grid_spec = GridSpec(dlat=0.05f0, dlon=0.05f0)
time_spec = TimeSpec(DateTime("2019-01-01"), DateTime("2019-12-31"), Dates.Day(1))

grid_center(config, grid_spec, time_spec;
            geo_provider=:modis,
            veg_indices=true, outfile="modis_2019.nc")
```

Experienced users can also use the dispatch-based library entry point:

```julia
grid(config, grid_spec, time_spec, SubpixelGridding(n_oversample=20);
     outfile="l2_output.nc")

grid(config, grid_spec, time_spec, CenterPointGridding();
     geo_provider=:modis, outfile="modis_output.nc")
```

To pre-generate one or more MODIS geolocation tiles:

```bash
julia --project=. bin/generate_modis_geolocation.jl \
    --tiles h08v04,h09v04 \
    --cacheDir ~/.cache/SatelliteGridding/modis/sinusoidal_2400px_v1
```

## Quick Start — Command Line

The CLI entry point is `bin/grid.jl`. It has two commands: **`l2`** (footprint oversampling) and **`center`** (MODIS-style).

```bash
# TROPOMI SIF — 8-day composites at 0.5°
julia --project=. bin/grid.jl l2 \
    --config examples/tropomi_sif.toml \
    --dLat 0.5 --dLon 0.5 --dDays 8 \
    --startDate 2019-01-01 --stopDate 2019-12-31 \
    -o tropomi_sif_2019.nc

# OCO-2 XCO₂ — monthly composites at 2°
julia --project=. bin/grid.jl l2 \
    --config examples/oco2_xco2.toml \
    --dLat 2.0 --dLon 2.0 --monthly --dDays 1 \
    --startDate 2015-01-01 --stopDate 2020-12-31 \
    -o oco2_xco2_monthly.nc

# TROPOMI NO₂ with KA CPU backend for parallel execution
julia --project=. bin/grid.jl l2 \
    --config examples/tropomi_no2.toml \
    --backend cpu --dLat 0.25 --dLon 0.25 \
    --startDate 2019-01-01 --stopDate 2019-12-31 \
    -o tropomi_no2.nc

# MODIS reflectance — center-coordinate mode with vegetation indices
julia --project=. bin/grid.jl center \
    --config examples/modis_reflectance.toml \
    --dLat 0.05 --dLon 0.05 --geoProvider modis --vegIndices \
    --startDate 2019-01-01 --stopDate 2019-12-31 \
    -o modis_2019.nc
```

Run `julia --project=. bin/grid.jl l2 --help` for all options.

### CLI flags (`l2` command)

| Flag | Default | Description |
|:-----|:--------|:------------|
| `--config`, `-c` | *required* | TOML (or JSON) configuration file |
| `--outFile`, `-o` | `gridded_output.nc` | Output NetCDF path |
| `--dLat` / `--dLon` | 1.0 | Grid resolution in degrees |
| `--latMin` / `--latMax` | -90 / 90 | Latitude bounds |
| `--lonMin` / `--lonMax` | -180 / 180 | Longitude bounds |
| `--startDate` / `--stopDate` | `2018-03-07` / `2018-10-31` | Date range (YYYY-MM-DD) |
| `--dDays` | 8 | Time step in days (or months with `--monthly`) |
| `--monthly` | false | Use months instead of days |
| `--nOversample` | 0 (auto) | Sub-pixel factor (0 = auto-compute) |
| `--compSTD` | false | Compute standard deviation |
| `--backend` | `sequential` | `sequential`, `cpu`, or `cuda` |

## Configuration Files

Data sources are defined in TOML files. Each config tells the package where to find input files, which NetCDF variables to grid, and what quality filters to apply:

```toml
filePattern = "S5P_PAL__L2B_SIF____YYYYMMDD*.nc"
folder = "/path/to/tropomi/data/"

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
```

### Sections

| Section | Required | Description |
|:--------|:---------|:------------|
| `[basic]` | required for `l2`; optional for `center` | Variable paths for `lat`, `lon`, `lat_bnd`, `lon_bnd`. Center mode only needs `lat`/`lon`, and MODIS can use generated sinusoidal geolocation instead |
| `[grid]` | yes | Output name → input variable path. All listed variables are gridded |
| `[filter]` | no | Quality filters using string expressions (see below) |
| `[center]` | no | Scale/fill options, MODIS pixel count, and explicit vegetation-index band roles for center gridding |
| `filePattern` | yes | Glob pattern with `YYYY`/`MM`/`DD`/`DOY`/`YYMMDD` placeholders |
| `folder` | yes | Root directory for input files (also supports date placeholders) |

### Filter expressions

Filters are defined in the `[filter]` section using intuitive string syntax:

| Expression | Meaning | Example |
|:-----------|:--------|:--------|
| `"< value"` | Less than | `"< 80"` |
| `"> value"` | Greater than | `"> 0.5"` |
| `"== value"` | Equality | `"== 0"` |
| `"lo < x < hi"` | Range (exclusive) | `"1600 < x < 2200"` |

All filters are combined with AND logic.

### Available example configs

| File | Instrument | Variables |
|:-----|:-----------|:----------|
| `examples/tropomi_sif.toml` | TROPOMI | SIF 743nm, SIF 735nm, cloud fraction |
| `examples/tropomi_no2.toml` | TROPOMI | Tropospheric NO₂ column |
| `examples/tropomi_ch4.toml` | TROPOMI | CH₄ mixing ratio |
| `examples/tropomi_co.toml` | TROPOMI | CO total column |
| `examples/oco2_sif.toml` | OCO-2 | SIF 757nm, SIF 771nm |
| `examples/oco2_xco2.toml` | OCO-2 | XCO₂, AOD, albedo |
| `examples/modis_reflectance.toml` | MODIS | Surface reflectance bands |

## Backends

The package supports three compute backends:

| Backend | CLI flag | Julia API | Description |
|:--------|:---------|:----------|:------------|
| Sequential | `--backend sequential` | default (`backend=nothing`) | Welford's online algorithm. Single-threaded. Supports `--compSTD`. |
| KA CPU | `--backend cpu` | `backend=CPU()` | KernelAbstractions CPU kernels. Parallel corner sorting + sub-pixel computation. Sum-based accumulation. |
| KA CUDA | `--backend cuda` | `backend=CUDABackend()` | Full GPU pipeline. Requires CUDA.jl (`using CUDA`). |

The KA backends use sum-based accumulation (`mean = sum/weight`) instead of Welford's incremental mean, which is fully parallelizable but does not support `--compSTD` in a single pass.

### GPU example

```julia
using CUDA, KernelAbstractions, SatelliteGridding, Dates

config = load_config("examples/tropomi_sif.toml")
grid_spec = GridSpec(dlat=0.5f0, dlon=0.5f0)
time_spec = TimeSpec(DateTime("2019-07-01"), DateTime("2019-07-31"), Dates.Day(16))

grid_l2(config, grid_spec, time_spec;
        backend=CUDABackend(), outfile="output.nc")
```

## How oversampling works

Each satellite footprint is a quadrilateral defined by 4 corner lat/lon pairs:

1. Corners are sorted into counter-clockwise (CCW) order
2. Two opposite edges are subdivided into n equally-spaced points
3. Corresponding points on opposite edges are connected and subdivided into n points
4. Result: n×n sub-pixel positions filling the footprint interior
5. Each sub-pixel is mapped to a grid cell via `floor()` and contributes weight `1/n²`

**Fast path**: When all 4 corners fall in the same grid cell, the full weight goes to that cell (no sub-pixel computation).

**Adaptive n**: By default, n is auto-computed from the footprint-to-grid-cell size ratio (~3 sub-pixels per grid cell width, clamped to [2, 20]). Override with `n_oversample=10` or `--nOversample 10`.

**Corner ordering**: Different satellite products store corners in different orders. The package normalizes them to CCW order automatically using a 5-comparator sorting network.

## Running tests

```julia
using Pkg
Pkg.test("SatelliteGridding")
```

Tests use synthetic data only — no real satellite files needed. Currently 220 tests covering oversampling, filtering, corner sorting, KA kernels, and I/O.

## Documentation

Full documentation is available at **[cfranken.github.io/SatelliteGridding](https://cfranken.github.io/SatelliteGridding/dev/)**.

Docs are built automatically via GitHub Actions on every push to `master`. To build locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated HTML will be in `docs/build/`.

## Legacy scripts

The original standalone scripts are preserved in `legacy/` for reference. The canonical script was `gridL2_Dates.jl`. Legacy JSON configuration files are in `legacy/jsonFiles/`.
