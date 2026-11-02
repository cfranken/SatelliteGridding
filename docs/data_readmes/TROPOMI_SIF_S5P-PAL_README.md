# TROPOMI SIF (S5P-PAL L2B) — Solar-induced chlorophyll fluorescence

> Drop this file at `/kiwi-data/Data/satellite/TROPOMI/TROPOMI_SIF_S5P-PAL/README.md`
> (needs write access to that user-owned directory; the existing short
> `README` already there can be removed/renamed).

## Dataset

Caltech retrieval of solar-induced chlorophyll fluorescence (SIF) at 743
and 735 nm from the Sentinel-5P/TROPOMI instrument, distributed through
the Sentinel-5P+ Innovation **Products Algorithm Laboratory (PAL)** data
portal. Daily Level-2B granules cover the full TROPOMI swath
(~2600 km wide, ~3.5 × 5 km footprints).

- Mission start: 2017-10 (commissioning); science data: **from 2018-05-01**
- Each NetCDF file is one orbital day (UTC), ~600 MB compressed
- Key variables (under `PRODUCT/`):
  - `SIF_743`, `SIF_735` — raw retrieved SIF (mW m⁻² sr⁻¹ nm⁻¹)
  - `SIF_Corr_743`, `SIF_Corr_735` — PAR/airmass-corrected SIF
  - `SIF_ERROR_743`, `SIF_ERROR_735` — single-sounding uncertainty
  - `latitude`, `longitude`, `latitude_bounds`, `longitude_bounds`
    (4-corner footprint polygons under `PRODUCT/SUPPORT_DATA/GEOLOCATIONS/`)
  - `cloud_fraction_L2`, `solar_zenith_angle`

## Source

| | |
|---|---|
| Portal | <https://data-portal.s5p-pal.com/products/troposif.html> |
| Algorithm | Köhler et al. (GRL 2018) Caltech retrieval |
| License | <https://sentinels.copernicus.eu/documents/247904/690755/Sentinel_Data_Legal_Notice> |
| Local maintainer | E. Metz (cron-refreshed daily ~05:00 UTC) |

## Citation

> Köhler, P., Frankenberg, C., Magney, T. S., Guanter, L., Joiner, J., &
> Landgraf, J. (2018). *Global retrievals of solar-induced chlorophyll
> fluorescence with TROPOMI: First results and intersensor comparison
> to OCO-2.* Geophysical Research Letters, 45(19), 10,456–10,463.
> <https://doi.org/10.1029/2018GL079031>

Sentinel-5P/TROPOMI mission paper:

> Veefkind, J. P., et al. (2012). *TROPOMI on the ESA Sentinel-5
> Precursor: A GMES mission for global observations of the atmospheric
> composition for climate, air quality and ozone layer applications.*
> Remote Sensing of Environment, 120, 70–83.
> <https://doi.org/10.1016/j.rse.2011.09.027>

## Layout

```
TROPOMI_SIF_S5P-PAL/
├── original/             # raw daily L2B NetCDFs from PAL (current)
│   └── S5P_PAL__L2B_SIF____<YYYYMMDDThhmmss>_<...>.nc
├── regridded/            # group-derived gridded products
│   └── daily_1x1/        # 1° daily filtered (used by various projects)
├── scripts/              # downloaders + maintenance
│   ├── Download_AllAvailableData.py
│   └── cron_script_TROPOMI_SIF_S5P-PAL.sh
└── README.md
```

## Download

The initial backfill used `scripts/Download_AllAvailableData.py`
(scrapes the PAL portal). New data is pulled by the daily cron job
`scripts/cron_script_TROPOMI_SIF_S5P-PAL.sh` at ~05:00 UTC, which also
publishes a global preview to
<https://web.gps.caltech.edu/~cfranken/tropomi-dashboard/>.

To re-fetch a subset manually (PAL allows unauthenticated bulk download):

```bash
DEST=/kiwi-data/Data/satellite/TROPOMI/TROPOMI_SIF_S5P-PAL/original
mkdir -p "$DEST"
# example: re-fetch 2024-08-01
DAY=20240801
curl -sS "https://data-portal.s5p-pal.com/api/s5p-l2b-sif/files?date=$DAY" \
  | jq -r '.[].downloadUrl' \
  | xargs -n1 -P4 wget -nc -P "$DEST"
```

## Gridded derivatives produced in this group

- `~/code/gitHub/Gridding/SatelliteGridding/tropomi_sif_16d_1deg_africa_doy_2018_2025.nc`
  — 16-day, 1°, Africa-box, DOY-anchored
- `regridded/daily_1x1/TROPOMI_sif_20180501_20251231_1deg_daily_gt-2SIF743lt5_SIF743ERRORlt10.nc`
  — 1° daily filtered (used in the SMAP–SIF analyses)

See `~/code/gitHub/Gridding/SatelliteGridding/bin/run_tropomi_l2_16day.sh`
for the 16-day footprint-oversampled pipeline (`grid.jl l2`).
