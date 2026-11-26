"""Generate independent HEALPix reference values with healpy (nested ordering).

Writes three CSV files into the directory given as argv[1]:
  healpix_ang2pix_ref.csv : nside, lon_deg, lat_deg, pix_nest
  healpix_pix2ang_ref.csv : nside, pix_nest, lon_deg, lat_deg   (cell centres)
  healpix_corners_ref.csv : nside, pix_nest, k, lon_deg, lat_deg (4 vertices, healpy order N,W,S,E)
"""
import csv
import sys

import healpy as hp
import numpy as np

out = sys.argv[1]
rng = np.random.default_rng(20261125)
nsides = [1, 2, 4, 64, 1024, 8192]

# Points: random uniform on sphere + special latitudes/longitudes.
n_rand = 400
u = rng.uniform(-1, 1, n_rand)
lat_r = np.degrees(np.arcsin(u))
lon_r = rng.uniform(-180, 180, n_rand)
special_lat = [90.0, -90.0, 0.0, np.degrees(np.arcsin(2 / 3)), -np.degrees(np.arcsin(2 / 3)),
               89.999, -89.999, 41.8, 41.82, -41.8, -41.82, 60.0, -75.0]
special_lon = [-180.0, 180.0, 0.0, 45.0, 90.0, -90.0, 135.0, -135.0, 179.9999, -179.9999, 10.0]
lat_s, lon_s = np.meshgrid(special_lat, special_lon)
lat_all = np.concatenate([lat_r, lat_s.ravel()])
lon_all = np.concatenate([lon_r, lon_s.ravel()])

with open(f"{out}/healpix_ang2pix_ref.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["nside", "lon_deg", "lat_deg", "pix_nest"])
    for nside in nsides:
        pix = hp.ang2pix(nside, lon_all, lat_all, nest=True, lonlat=True)
        for lo, la, p in zip(lon_all, lat_all, pix):
            w.writerow([nside, repr(float(lo)), repr(float(la)), int(p)])

with open(f"{out}/healpix_pix2ang_ref.csv", "w", newline="") as f, \
     open(f"{out}/healpix_corners_ref.csv", "w", newline="") as g:
    w = csv.writer(f)
    wc = csv.writer(g)
    w.writerow(["nside", "pix_nest", "lon_deg", "lat_deg"])
    wc.writerow(["nside", "pix_nest", "k", "lon_deg", "lat_deg"])
    for nside in nsides:
        npix = hp.nside2npix(nside)
        if npix <= 200:
            pix = np.arange(npix)
        else:
            pix = np.unique(np.concatenate([
                rng.integers(0, npix, 150),
                [0, 1, 2, 3, npix - 1, npix // 2, npix // 12, 4 * npix // 12, 8 * npix // 12],
            ]))
        lon, lat = hp.pix2ang(nside, pix, nest=True, lonlat=True)
        for p, lo, la in zip(pix, lon, lat):
            w.writerow([nside, int(p), repr(float(lo)), repr(float(la))])
        vec = hp.boundaries(nside, pix, step=1, nest=True)  # (n, 3, 4)
        vec = np.atleast_3d(vec) if vec.ndim == 2 else vec
        for i, p in enumerate(pix):
            for k in range(4):
                x, y, z = vec[i, :, k]
                lo, la = hp.vec2ang(np.array([x, y, z]), lonlat=True)
                wc.writerow([nside, int(p), k + 1, repr(float(lo[0])), repr(float(la[0]))])

print("healpy", hp.__version__, "points", len(lat_all))
