# CLI Reference

`SatelliteGridding.jl` includes a command-line interface via `bin/grid.jl`.

## Usage

```bash
julia --project=. bin/grid.jl <command> [options]
```

## Commands

### `l2` — Grid Level-2 Data

Grid satellite Level-2 data with footprint oversampling.

```bash
julia --project=. bin/grid.jl l2 [options]
```

#### Options

| Flag | Type | Default | Description |
|:-----|:-----|:--------|:------------|
| `--config`, `-c` | String | *required* | TOML/JSON configuration file |
| `--outFile`, `-o` | String | `gridded_output.nc` | Output NetCDF filename |
| `--latMin` | Float32 | -90 | Lower latitude bound |
| `--latMax` | Float32 | 90 | Upper latitude bound |
| `--lonMin` | Float32 | -180 | Lower longitude bound |
| `--lonMax` | Float32 | 180 | Upper longitude bound |
| `--dLat` | Float32 | 1.0 | Latitude resolution (degrees) |
| `--dLon` | Float32 | 1.0 | Longitude resolution (degrees) |
| `--startDate` | String | `2018-03-07` | Start date (YYYY-MM-DD) |
| `--stopDate` | String | `2018-10-31` | Stop date (YYYY-MM-DD) |
| `--dDays` | Int | 8 | Time step in days (or months with `--monthly`) |
| `--monthly` | Flag | false | Use months instead of days for time step |
| `--oversample_temporal` | Float32 | 1.0 | Temporal oversampling factor |
| `--nOversample` | Int | 0 (auto) | Sub-pixel factor (0 = auto-compute) |
| `--compSTD` | Flag | false | Compute standard deviation |
| `--backend` | String | `sequential` | Compute backend: `sequential`, `cpu`, or `cuda` |
| `--keepGoing` | Flag | false | Continue after per-file processing errors |

#### Examples

Global TROPOMI SIF at 0.5° / 8-day composites:
```bash
julia --project=. bin/grid.jl l2 \
    --config examples/tropomi_sif.toml \
    --dLat 0.5 --dLon 0.5 \
    --startDate 2019-01-01 --stopDate 2019-12-31 \
    --dDays 8 \
    -o tropomi_sif_2019_8day_05deg.nc
```

Regional OCO-2 XCO₂ with KA CPU backend:
```bash
julia --project=. bin/grid.jl l2 \
    --config examples/oco2_xco2.toml \
    --latMin -60 --latMax 80 \
    --dLat 2.0 --dLon 2.0 \
    --startDate 2019-01-01 --stopDate 2019-12-31 \
    --dDays 30 \
    --backend cpu \
    -o oco2_xco2_2019.nc
```

Monthly composites with standard deviation:
```bash
julia --project=. bin/grid.jl l2 \
    --config examples/tropomi_no2.toml \
    --dLat 0.25 --dLon 0.25 \
    --startDate 2019-01-01 --stopDate 2019-12-31 \
    --monthly --dDays 1 \
    --compSTD \
    -o tropomi_no2_2019_monthly.nc
```

### `center` — Grid Center-Coordinate Data

Grid data using center coordinates only (no footprint bounds). Suitable for
MODIS-style data where each pixel maps to exactly one grid cell.

```bash
julia --project=. bin/grid.jl center [options]
```

#### Options

| Flag | Type | Default | Description |
|:-----|:-----|:--------|:------------|
| `--config`, `-c` | String | *required* | TOML/JSON configuration file |
| `--outFile`, `-o` | String | `gridded_output.nc` | Output NetCDF filename |
| `--latMin` | Float32 | -90 | Lower latitude bound |
| `--latMax` | Float32 | 90 | Upper latitude bound |
| `--lonMin` | Float32 | -180 | Lower longitude bound |
| `--lonMax` | Float32 | 180 | Upper longitude bound |
| `--dLat` | Float32 | 0.5 | Latitude resolution (degrees) |
| `--dLon` | Float32 | 0.5 | Longitude resolution (degrees) |
| `--startDate` | String | `2018-01-01` | Start date (YYYY-MM-DD) |
| `--stopDate` | String | `2018-12-31` | Stop date (YYYY-MM-DD) |
| `--dDays` | Int | 1 | Time step in days (or months with `--monthly`) |
| `--monthly` | Flag | false | Use months instead of days for time step |
| `--geoProvider` | String | `auto` | Center geolocation provider: `auto`, `variables`, `lut`, or `modis` |
| `--geoTable` | String | *(none)* | Path to legacy monolithic geolocation lookup table (NetCDF) |
| `--geoCache` | String | user cache | Directory for generated per-tile MODIS sinusoidal geolocation |
| `--vegIndices` | Flag | false | Compute vegetation indices (EVI, NDVI, NIRv, NDWI) |
| `--keepGoing` | Flag | false | Continue after per-file processing errors |

#### Example

```bash
julia --project=. bin/grid.jl center \
    --config examples/modis_reflectance.toml \
    --dLat 0.05 --dLon 0.05 \
    --startDate 2019-01-01 --stopDate 2019-12-31 \
    --geoProvider modis \
    --vegIndices \
    -o modis_2019.nc
```

MODIS geolocation tiles are generated on demand into the user cache. They can
also be generated ahead of time:

```bash
julia --project=. bin/generate_modis_geolocation.jl \
    --tiles h08v04,h09v04 \
    --cacheDir ~/.cache/SatelliteGridding/modis/sinusoidal_2400px_v1
```

## Backends

| Backend | Flag | Description |
|:--------|:-----|:------------|
| Sequential | `--backend sequential` | Default. Uses Welford's online algorithm for mean/std. Single-threaded. Supports `--compSTD`. |
| KA CPU | `--backend cpu` | KernelAbstractions CPU backend. Parallel sort + subpixel computation. Sum-based accumulation. |
| KA CUDA | `--backend cuda` | KernelAbstractions GPU backend. Requires CUDA.jl. All computation on GPU with atomic scatter. |

The `cpu` and `cuda` backends use sum-based accumulation (mean = sum/weight) instead of
Welford's incremental mean. This is fully parallelizable but does not support `--compSTD`
in a single pass. Standard deviation with KA backends requires a two-pass approach.
