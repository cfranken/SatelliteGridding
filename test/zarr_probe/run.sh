#!/usr/bin/env bash
#
# Mandatory ordering verification (spec §4): write a ramp probe through the real
# Julia writer, then assert the index→coordinate mapping in BOTH xarray (numpy
# consumers) and zarrita (the browser dashboard). Any transpose / lat-flip /
# lon-roll / C-vs-F mistake fails loudly here.
#
# Prereqs (one-time):
#   python3 -m pip install --user zarr        # xarray.open_zarr needs zarr
#   (cd test/zarr_probe && npm install zarrita numcodecs)
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
JL=${JULIA:-/home/cfranken/.julia/juliaup/julia-1.12.5+0.x64.linux.gnu/bin/julia}
STORE=probe.zarr
PORT=${PORT:-8123}

rm -rf "$STORE"
echo "[1/4] writing ramp probe via the real Julia writer ..."
"$JL" --project="$ROOT" probe.jl "$STORE"

echo "[2/4] Python / xarray ordering check ..."
python3 check.py "$STORE"

echo "[3/4] serving store over HTTP for zarrita ..."
python3 -m http.server "$PORT" >/dev/null 2>&1 &
SRV=$!
trap 'kill "$SRV" 2>/dev/null || true; rm -rf "'"$STORE"'"' EXIT
sleep 1

echo "[4/4] Node / zarrita ordering check ..."
node check.mjs "http://localhost:$PORT/$STORE"

echo "ALL ORDERING CHECKS PASSED"
