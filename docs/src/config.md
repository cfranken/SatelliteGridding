# Configuration Reference

Data source configurations are defined in TOML files (legacy JSON is also supported).
Each config file tells `SatelliteGridding.jl` where to find input data, which variables
to grid, and what quality filters to apply.

## File Structure

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
```

## Top-Level Fields

| Field | Description |
|:------|:------------|
| `filePattern` | Glob pattern for input files. Supports date placeholders: `YYYY` (year), `MM` (month), `DD` (day), `DOY` (day of year) |
| `folder` | Root directory for input data. Also supports `YYYY`/`MM`/`DD` placeholders for date-organized directory trees |

### Date Placeholders

The file pattern and folder path support these placeholders:

- `YYYY` — 4-digit year (e.g., `2019`)
- `YYMMDD` — compact date (e.g., `190701`)
- `MM` — 2-digit month (e.g., `07`)
- `DD` — 2-digit day (e.g., `01`)
- `DOY` — 3-digit day of year (e.g., `182`)

Example for OCO-2 with date-organized directories:
```toml
filePattern = "oco2_LtCO2_YYMMDD_*.nc4"
folder = "/data/OCO2/B9_SIF/"
```

## `[basic]` Section (Required)

Maps internal coordinate keys to NetCDF variable paths.

| Key | Description |
|:----|:------------|
| `lat` | Center latitude variable |
| `lon` | Center longitude variable |
| `lat_bnd` | Corner latitude bounds (N×4 or 4×N array) |
| `lon_bnd` | Corner longitude bounds (N×4 or 4×N array) |

For center-coordinate gridding (`grid_center`), only `lat` and `lon` are needed.

NetCDF group paths use `/` separators (e.g., `"PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds"`).

## `[grid]` Section (Required)

Maps output variable names to input NetCDF variable paths. The key becomes the
variable name in the output NetCDF file.

```toml
[grid]
sif_743 = "PRODUCT/SIF_743"           # Output name = "sif_743"
cloud_fraction = "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction_L2"
```

Variables are gridded in the order they appear. For MODIS data, band order matters
for vegetation index computation.

## `[filter]` Section (Optional)

Quality filter rules using intuitive string expressions. Each key is a NetCDF
variable path, and the value is a filter expression.

### Filter Expressions

| Syntax | Meaning | Example |
|:-------|:--------|:--------|
| `"< value"` | Less than | `"< 80"` |
| `"> value"` | Greater than | `"> 0.5"` |
| `"== value"` | Equality | `"== 0"` |
| `"lo < x < hi"` | Range (exclusive) | `"1600 < x < 2200"` |

All filters are combined with AND logic — an observation must pass all filters
to be included in the gridding.

### Examples

```toml
[filter]
# Solar zenith angle less than 80 degrees
"PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle" = "< 80"

# Quality flag equals 0 (good quality)
"xco2_quality_flag" = "== 0"

# Methane in valid range
"PRODUCT/methane_mixing_ratio_bias_corrected" = "1600 < x < 2200"
```

## Example Configurations

### TROPOMI SIF

```toml
filePattern = "S5P_PAL__L2B_SIF____YYYYMMDD*.nc"
folder = "/data/TROPOMI/SIF/"

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
```

### OCO-2 XCO₂

```toml
filePattern = "oco2_LtCO2_YYMMDD_*.nc4"
folder = "/data/OCO2/B9_SIF/"

[basic]
lat = "latitude"
lon = "longitude"
lat_bnd = "vertex_latitude"
lon_bnd = "vertex_longitude"

[grid]
xco2 = "xco2"
aod = "Retrieval/aod_total"

[filter]
"xco2_quality_flag" = "== 0"
```

### TROPOMI CH₄

```toml
filePattern = "S5P_RPRO_L2__CH4____YYYYMMDD*.nc"
folder = "/data/TROPOMI/CH4/"

[basic]
lat = "PRODUCT/latitude"
lon = "PRODUCT/longitude"
lat_bnd = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds"
lon_bnd = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds"

[grid]
methane = "PRODUCT/methane_mixing_ratio_bias_corrected"

[filter]
"PRODUCT/qa_value" = "> 0.5"
"PRODUCT/methane_mixing_ratio_bias_corrected" = "1600 < x < 2200"
```

## Legacy JSON Format

The original JSON format is still supported:

```json
{
    "basic": {
        "lat": "PRODUCT/latitude",
        "lon": "PRODUCT/longitude",
        "lat_bnd": "...",
        "lon_bnd": "..."
    },
    "grid": {
        "sif_743": "PRODUCT/SIF_743"
    },
    "filter_lt": {
        "PRODUCT/.../solar_zenith_angle": 80
    },
    "filter_gt": {
        "PRODUCT/qa_value": 0.5
    },
    "filePattern": "...",
    "folder": "..."
}
```
