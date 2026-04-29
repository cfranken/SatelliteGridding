# # Tutorial: Gridding TROPOMI SIF Data
#
# This tutorial walks through gridding TROPOMI Solar-Induced Fluorescence (SIF)
# data onto a regular lat/lon grid using `SatelliteGridding.jl`.
#
# ## Setup
#
# First, load the package and define the time range:

using SatelliteGridding
using Dates

# ## Step 1: Load the Configuration
#
# Configuration files (TOML format) tell the package where to find data files,
# which NetCDF variables to read, and what quality filters to apply.
#
# Here's what a typical TROPOMI SIF config looks like:
#
# ```toml
# filePattern = "S5P_PAL__L2B_SIF____YYYYMMDD*.nc"
# folder = "/path/to/tropomi/sif/"
#
# [basic]
# lat = "PRODUCT/latitude"
# lon = "PRODUCT/longitude"
# lat_bnd = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds"
# lon_bnd = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds"
#
# [grid]
# sif_743 = "PRODUCT/SIF_743"
# sif_735 = "PRODUCT/SIF_735"
#
# [filter]
# "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle" = "< 80"
# ```

config = load_config(joinpath(pkgdir(SatelliteGridding), "examples", "tropomi_sif.toml"))

# The `config` object contains:
# - `basic`: Coordinate variable paths (lat, lon, and corner bounds)
# - `grid_vars`: Variables to be gridded (e.g., SIF at 743nm and 735nm)
# - `filters`: Quality filter rules (here: solar zenith angle < 80°)
# - `file_pattern`: Glob pattern with date placeholders
# - `folder`: Root data directory

# ## Step 2: Define the Output Grid
#
# A `GridSpec` defines the spatial extent and resolution of the output grid.
# All values are in degrees. The package uses Float32 for memory efficiency.

grid_spec = GridSpec(
    lat_min = -60f0,    # Southern bound
    lat_max =  80f0,    # Northern bound
    lon_min = -180f0,   # Western bound
    lon_max =  180f0,   # Eastern bound
    dlat = 0.5f0,       # Latitude resolution (degrees)
    dlon = 0.5f0        # Longitude resolution (degrees)
)

# This creates a 720×280 grid. The `grid_spec.lat` and `grid_spec.lon` vectors
# contain cell center coordinates.
println("Grid dimensions: $(length(grid_spec.lon)) × $(length(grid_spec.lat))")

# ## Step 3: Define the Time Specification
#
# A `TimeSpec` controls temporal binning. You can use `Dates.Day` or `Dates.Month`
# for the time step.

time_spec = TimeSpec(
    DateTime("2019-07-01"),     # Start date
    DateTime("2019-07-31"),     # Stop date
    Dates.Day(16)               # 16-day composites
)

# For monthly composites:
# ```julia
# time_spec = TimeSpec(DateTime("2019-01-01"), DateTime("2019-12-31"),
#                      Dates.Month(1))
# ```

# ## Step 4: Run the Gridding
#
# The `grid_l2` function processes all matching files, applies filters,
# subdivides footprints, and writes the output NetCDF file.

## grid_l2(config, grid_spec, time_spec;
##         outfile = "tropomi_sif_july2019.nc",
##         n_oversample = 10)    # 10×10 sub-pixels per footprint

# ### Key Options
#
# - `n_oversample`: Sub-pixel subdivision factor. Higher = more accurate spatial
#   distribution but slower. Default: auto-computed from footprint/grid ratio.
# - `compute_std`: Set to `true` to also compute per-cell standard deviation.
# - `backend`: Compute backend — `nothing` (sequential Welford), `CPU()` (KA parallel),
#   or `CUDABackend()` (GPU).
#
# ### Using the KA Backend
#
# For large datasets, the KernelAbstractions backend parallelizes the computation:
#
# ```julia
# using KernelAbstractions
# grid_l2(config, grid_spec, time_spec;
#         outfile = "output.nc",
#         backend = CPU())
# ```
#
# The KA backend uses sum-based accumulation (instead of Welford's incremental mean),
# which is fully parallelizable. The mean is computed at the end: `mean = sum / weight`.

# ## Step 5: Inspect the Output
#
# The output NetCDF file contains:
# - `lat`, `lon`: Grid cell center coordinates
# - `time`: Time step dates
# - `n`: Number of observations per cell
# - One variable per entry in `[grid]` (e.g., `sif_743`, `sif_735`)
# - Optional `_std` suffix variables if `compute_std=true`
#
# ```julia
# using NCDatasets
# ds = Dataset("tropomi_sif_july2019.nc")
# sif = ds["sif_743"][:, :, 1]   # First time step
# weights = ds["n"][:, :, 1]
# close(ds)
# ```

# ## Center-Coordinate Gridding (MODIS)
#
# For instruments without footprint bounds (e.g., MODIS), use `grid_center`:
#
# ```julia
# config = load_config("examples/modis_reflectance.toml")
# grid_spec = GridSpec(dlat=0.05f0, dlon=0.05f0)
# time_spec = TimeSpec(DateTime("2019-01-01"), DateTime("2019-12-31"),
#                      Dates.Day(1))
# grid_center(config, grid_spec, time_spec;
#             geo_provider=:modis,
#             veg_indices=true,
#             outfile="modis_2019.nc")
# ```
