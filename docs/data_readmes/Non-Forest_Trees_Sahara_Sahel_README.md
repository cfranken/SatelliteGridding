# Sahel/Sahara tree carbon density — 100 m semi-arid tree mapping

> Drop this file at
> `/kiwi-data/Data/satellite/Non-Forest_Trees_Sahara_Sahel/README.md`.
> The directory is owned by egreich — owner needs to install or `chmod`.

## Dataset

A near-continuous 100 m raster of **above-ground tree carbon density**
across the West African Sahel and Sahara, derived by Tucker, Brandt and
colleagues from sub-meter Maxar satellite imagery (2017–2020 epoch).
Each ~50 cm pixel within the original tree-crown classification was
allometrically converted to a per-crown carbon mass, then aggregated to
a 100 m grid of "tree carbon density" (Mg C ha⁻¹).

- Static (single multi-year epoch, not a time series)
- Coverage: lon ~−18° to +42°, lat ~9.4° to +24.1° (Sahel/Sahara band only)
- 5 raster bands in the GeoTIFF (band 1 = tree carbon density; later
  bands contain uncertainty and intermediate products documented in
  `docs/`)
- Format: cloud-optimized GeoTIFF, EPSG:4326-like geographic CRS,
  float32 with nodata = −1

## Source

| | |
|---|---|
| Distribution portal | ORNL DAAC, dataset ID 2117 |
| Identifier | `Semi-Arid_Tree_Carbon_50cm_2117_1` |
| Data viewer description | `docs/Data_viewer_description_Carbon_Stocks_of_Indiv_Trees.pdf` |
| Algorithm summary | `docs/Semi-Arid_Tree_Carbon_50cm.pdf` |

## Citation

> Tucker, C. J., Brandt, M., Hiernaux, P., Kariryaa, A., Rasmussen, K.,
> Small, J., Igel, C., Reiner, F., Melocik, K., Meyer, J., Sinno, S.,
> Romero, E., Glennie, E., Fitts, Y., Morin, A., Pinzon, J., McClain, D.,
> Rubin, P. R., Manchego, L. C., Carr, J., Dziubek, M., & Fensholt, R.
> (2023). *Sub-continental-scale carbon stocks of individual trees in
> African drylands.* Nature, 615(7950), 80–86.
> <https://doi.org/10.1038/s41586-022-05653-6>

Original tree-detection method paper:

> Brandt, M., et al. (2020). *An unexpectedly large count of trees in
> the West African Sahara and Sahel.* Nature, 587(7832), 78–82.
> <https://doi.org/10.1038/s41586-020-2824-5>

## Layout

```
Non-Forest_Trees_Sahara_Sahel/
├── data/
│   └── Semi-Arid_Tree_Carbon_50cm_2117_1-<timestamp>/
│       ├── wgs_84_mosaic_100m_tree_carbon_density.tif    # 8.7 GB
│       └── wgs_84_mosaic_100m_tree_carbon_density.tif.sha256
└── docs/
    ├── Data_viewer_description_Carbon_Stocks_of_Indiv_Trees.pdf
    └── Semi-Arid_Tree_Carbon_50cm.pdf
```

## Download

Dataset is hosted at ORNL DAAC:

```bash
# requires Earthdata Login (~/.netrc with urs.earthdata.nasa.gov)
DEST=/kiwi-data/Data/satellite/Non-Forest_Trees_Sahara_Sahel/data
mkdir -p "$DEST"
# replace <SUBPATH> with the current ORNL listing for dataset 2117
wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies \
     --keep-session-cookies --no-check-certificate --auth-no-challenge=on \
     -P "$DEST" \
     "https://daac.ornl.gov/daacdata/agroclimatology/Sub-continental-scale_Carbon_stocks/<SUBPATH>/wgs_84_mosaic_100m_tree_carbon_density.tif"
```

Or via `earthaccess`:

```python
import earthaccess
earthaccess.login(strategy="netrc")
gs = earthaccess.search_data(
    short_name="Sub-continental-scale_Carbon_stocks_of_Individual_Trees_in_African_Drylands_2117")
earthaccess.download(gs, "/kiwi-data/Data/satellite/Non-Forest_Trees_Sahara_Sahel/data/")
```

## Reading and aggregating

The TIF is 66499 × 16145 pixels at 100 m, **8.7 GB on disk**. Use
windowed reads (`rasterio.windows.from_bounds`) when subsetting; full
loads are wasteful.

A 1° aggregator that resamples band 1 to a chosen target grid lives in
`~/code/gitHub/Gridding/SatelliteGridding/bin/grid_tree_density_africa.py`
— output is `tree_carbon_density_1deg_africa.nc` (alongside the gridded
SIF/SMAP/MODIS files in that repo). The aggregator uses simple per-cell
mean of the 100 m pixels — appropriate because the 100 m grid is
already a density (Mg C ha⁻¹).
