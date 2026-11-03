# Installing the staged README files

Three of the per-dataset READMEs need to be placed in directories that
are owned by other users on the kiwi host. To install them, the owner of
the directory (or root) needs to run:

```bash
# replace OWNER with the dir owner where needed
sudo install -o "$OWNER" -g cfranken -m 0644 \
    /home/cfranken/code/gitHub/Gridding/SatelliteGridding/docs/data_readmes/TROPOMI_SIF_S5P-PAL_README.md \
    /kiwi-data/Data/satellite/TROPOMI/TROPOMI_SIF_S5P-PAL/README.md

sudo install -o "$OWNER" -g cfranken -m 0644 \
    /home/cfranken/code/gitHub/Gridding/SatelliteGridding/docs/data_readmes/SMAP_README.md \
    /kiwi-data/Data/satellite/SMAP/README.md

sudo install -o "$OWNER" -g cfranken -m 0644 \
    /home/cfranken/code/gitHub/Gridding/SatelliteGridding/docs/data_readmes/Non-Forest_Trees_Sahara_Sahel_README.md \
    /kiwi-data/Data/satellite/Non-Forest_Trees_Sahara_Sahel/README.md
```

Or — if the directory owner is doing it themselves — just `cp`:

```bash
cp /home/cfranken/code/gitHub/Gridding/SatelliteGridding/docs/data_readmes/TROPOMI_SIF_S5P-PAL_README.md \
   /kiwi-data/Data/satellite/TROPOMI/TROPOMI_SIF_S5P-PAL/README.md
# (etc.)
```

The four READMEs that were already writable have been placed in-tree:

- `/kiwi-data/Data/README.md`
- `/kiwi-data/Data/satellite/MODIS/MCD43A4/README.md`
- `/kiwi-data/Data/satellite/GRACE/README.md`
- `/kiwi-data/Data/satellite/GRACE/JPL_RL06.3_MASCON_CRI/README.md`
- `/kiwi-data/Data/climate_indices/README.md`
