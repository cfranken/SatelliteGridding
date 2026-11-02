# SMAP L3 Enhanced Passive Soil Moisture (SPL3SMP_E v6)

> Drop this file at `/kiwi-data/Data/satellite/SMAP/README.md`
> (overwrites the shorter existing one). The directory is owned by
> egreich — use `chmod g+w` first or have the owner replace.

## Dataset

NASA SMAP (Soil Moisture Active Passive) mission Level-3 daily composite
of the **enhanced** L-band radiometer soil moisture retrieval, on the
9 km global EASE-Grid 2.0. One HDF5 file per day, separate AM (descending
overpass, ~06:00 LST) and PM (ascending, ~18:00 LST) groups inside each
file.

- Mission start: SMAP launched 2015-01-31; L3 enhanced product begins
  **2015-03-31**.
- Spatial resolution: 9 km EASE-Grid; "Enhanced" means the SCA-V
  brightness-temperature-based retrieval is posted on the 9 km grid
  using the Backus–Gilbert oversampling of the native 36 km footprint.
- Variable for analysis: `Soil_Moisture_Retrieval_Data_AM/soil_moisture`
  (m³/m³). Use AM only (literature consensus: AM correlates better with
  surface SM than PM due to thermal stratification effects).

## Source

| | |
|---|---|
| Product short name | `SPL3SMP_E` (v006) |
| Collection ID | `SPL3SMP_E.006` |
| DOI | <https://doi.org/10.5067/4DQ54OUIJ9DL> |
| NSIDC landing page | <https://nsidc.org/data/spl3smp_e/versions/6> |

## Citation

> O'Neill, P. E., Chan, S., Njoku, E. G., Jackson, T., Bindlish, R.,
> Chaubell, J., & Colliander, A. (2024). *SMAP L3 Radiometer Global Daily
> 9 km EASE-Grid Soil Moisture, Version 6.* NASA National Snow and Ice
> Data Center Distributed Active Archive Center, Boulder, Colorado USA.
> doi:10.5067/4DQ54OUIJ9DL.

Mission overview:

> Entekhabi, D., et al. (2010). *The Soil Moisture Active Passive (SMAP)
> mission.* Proceedings of the IEEE, 98(5), 704–716.
> <https://doi.org/10.1109/JPROC.2010.2043918>

## Layout

```
SMAP/
├── SMAP_L3_SM_P_E/                 # raw HDF5 archive (one file per day)
│   ├── SMAP_L3_SM_P_E_20150331_R19240_001.h5
│   ├── SMAP_L3_SM_P_E_20150331_R19240_001.h5.iso.xml
│   ├── SMAP_L3_SM_P_E_20150331_R19240_001.qa
│   └── ...
├── SMAP_L3_SM_P_E_netcdf/          # subset converted to NetCDF (partial)
└── nsidc-download_SPL3SMP_E.006_<date>.py    # NSIDC-generated puller
```

Year-by-year file count in the H5 archive:

| year | count |
|---|---|
| 2015 | 274 (from Mar 31) |
| 2016 | 364 |
| 2017 | 365 |
| 2018 | 365 |
| 2019 | 332 |
| 2020 | 366 |
| 2021 | 365 |
| 2022 | 319 |
| 2023 | 365 |
| 2024 | 366 |
| 2025 | 312 (through Nov 1) |

## Download

NSIDC distributes via the Earthdata DAAC; auth is your `~/.netrc`
machine-`urs.earthdata.nasa.gov` entry. Use either the NSIDC-generated
puller (`nsidc-download_SPL3SMP_E.006_*.py` in this dir) or
`earthaccess`:

```python
import earthaccess
earthaccess.login(strategy="netrc")
gs = earthaccess.search_data(
    short_name="SPL3SMP_E", version="006",
    temporal=("2025-01-01", "2025-12-31"))
earthaccess.download(gs, "/kiwi-data/Data/satellite/SMAP/SMAP_L3_SM_P_E/")
```

Or a single day with curl:

```bash
DAY=20240601
URL="https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP_E.006/$(date -d $DAY +%Y.%m.%d)/SMAP_L3_SM_P_E_${DAY}_R19240_001.h5"
curl -O -b ~/.urs_cookies -c ~/.urs_cookies -L -n "$URL"
```

## Reading

The HDF5 file has separate `Soil_Moisture_Retrieval_Data_AM` and
`..._PM` groups. The relevant fields per group:

| field | shape | description |
|---|---|---|
| `soil_moisture` | (latitude, longitude) | 0.02–0.5 m³/m³, fill = −9999 |
| `latitude`, `longitude` | 2D on EASE-Grid | per-pixel coordinates |
| `retrieval_qual_flag` | (latitude, longitude) | bit-packed QA |

## Gridded derivatives produced in this group

Group scripts re-bin SMAP AM to 1° / 16-day for the Africa analyses:

- `~/code/gitHub/Gridding/SatelliteGridding/smap_sm_16d_1deg_africa_doy_2018_2025.nc`
- `~/code/gitHub/Gridding/SatelliteGridding/smap_sm_16d_1deg_africa_doy_2015_2025.nc`

Gridder: `bin/grid_smap_africa_16day.py` (reads HDF5 directly; configurable
via `--start`/`--end`/`--ddays`).
