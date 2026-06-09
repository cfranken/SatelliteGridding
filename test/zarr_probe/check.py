#!/usr/bin/env python3
"""Ordering check via xarray/zarr (spec §4.2) — authoritative for numpy consumers."""
import sys

import numpy as np
import xarray as xr

path = sys.argv[1] if len(sys.argv) > 1 else "probe.zarr"
ds = xr.open_zarr(path, consolidated=True)

a = ds["probe"].isel(time=0).values            # must be (lat, lon) = (180, 360)
assert a.shape == (180, 360), f"shape {a.shape} != (180, 360)"
assert ds["probe"].dims == ("time", "lat", "lon"), ds["probe"].dims

for iy, ix in [(0, 0), (44, 59), (89, 180), (179, 359)]:
    got = a[iy, ix]
    exp = iy * 1000 + ix
    assert got == exp, f"@(iy={iy}, ix={ix}): {got} != {exp}"   # north-up, lon-east, no transpose

assert np.isclose(ds.lat.values[0], 89.5) and np.isclose(ds.lat.values[-1], -89.5), ds.lat.values[[0, -1]]
assert np.isclose(ds.lon.values[0], -179.5) and np.isclose(ds.lon.values[-1], 179.5), ds.lon.values[[0, -1]]

print("dims:", ds["probe"].dims, "shape:", ds["probe"].shape, "encoding:", ds["probe"].encoding.get("chunks"))
print("PYTHON ORDERING OK")
