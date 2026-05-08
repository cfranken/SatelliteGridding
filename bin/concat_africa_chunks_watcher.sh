#!/usr/bin/env bash
#
# One-shot watcher: waits for the Africa MODIS + TROPOMI Julia gridding
# children to exit, then runs ncrcat to concatenate each product's per-chunk
# netCDFs into the final output. Each product is handled independently, so
# whichever finishes first gets concatenated first.
#
# Designed to be launched detached (nohup ... &). Avoids `set -e` on purpose
# so a single concat failure doesn't take down the other product's watcher.

REPO="/home/cfranken/code/gitHub/Gridding/SatelliteGridding"
cd "$REPO" || exit 2

LOG="$REPO/africa_concat_watcher.log"
exec > "$LOG" 2>&1

echo "[$(date '+%Y-%m-%d %H:%M:%S')] watcher start, PID=$$"

watch_and_concat() {
    local product="$1"
    local pattern="$2"
    local outfile="$3"
    local chunks_dir="$4"

    while pgrep -f "$pattern" > /dev/null 2>&1; do
        sleep 60
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $product: all children exited; concatenating"

    local chunks=( "$chunks_dir"/chunk_*.nc )
    # Bash returns the literal glob if no match, so guard against that.
    if (( ${#chunks[@]} == 0 )) || [[ ! -f "${chunks[0]}" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $product: no chunk files found in $chunks_dir, skipping"
        return
    fi

    if ncrcat -O -h "${chunks[@]}" "$outfile"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $product OK -> $outfile  ($(du -h "$outfile" | cut -f1))"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $product FAILED ncrcat"
    fi
}

watch_and_concat "MODIS" \
    "grid.jl center.*africa_2018_2025" \
    "modis_mcd43a4_16d_1deg_africa_2018_2025.nc" \
    "modis_mcd43a4_16d_1deg_africa_2018_2025_chunks" &

watch_and_concat "TROPOMI" \
    "grid.jl l2.*tropomi_sif_16d_1deg_africa" \
    "tropomi_sif_16d_1deg_africa_2018_2025.nc" \
    "tropomi_sif_16d_1deg_africa_2018_2025_chunks" &

wait
echo "[$(date '+%Y-%m-%d %H:%M:%S')] watcher done"
