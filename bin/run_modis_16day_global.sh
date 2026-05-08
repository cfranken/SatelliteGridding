#!/usr/bin/env bash
#
# Grid MCD43A4 surface reflectance to a 1°, 16-day global product, 2018-2025.
# Outputs Nadir reflectance bands 1-7 plus EVI, NDVI, NIRv, NDWI vegetation
# indices (NDWI is co-produced; only NDVI/NIRv/EVI were requested).
#
# Fans out across time chunks as independent Julia processes (the center
# gridder is single-threaded), then concatenates the per-chunk netCDFs with
# ncrcat. Chunks are aligned to 16-day boundaries from START_DATE, so the
# concatenated output has the exact same cadence as a single-process run.
#
# Run from the SatelliteGridding repo root:
#   bash bin/run_modis_16day_global.sh
#
# Tunables (env vars):
#   JOBS=8                  number of concurrent Julia processes
#   OUTFILE=...             final concatenated netCDF
#   CHUNK_DIR=...           directory for per-chunk intermediates
#   START_DATE=2018-01-01
#   STOP_DATE=2025-12-31
#   DDAYS=16
#   DLAT=1.0   DLON=1.0
#   GEO_CACHE=...           pre-populated MODIS sinusoidal cache
#   KEEP_CHUNKS=1           keep per-chunk files after concat (default: delete)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

CONFIG="${CONFIG:-examples/modis_reflectance_kiwi.toml}"
OUTFILE="${OUTFILE:-modis_mcd43a4_16d_1deg_2018_2025.nc}"
START_DATE="${START_DATE:-2018-01-01}"
STOP_DATE="${STOP_DATE:-2025-12-31}"
DDAYS="${DDAYS:-16}"
DLAT="${DLAT:-1.0}"
DLON="${DLON:-1.0}"
JOBS="${JOBS:-8}"
GEO_CACHE="${GEO_CACHE:-$HOME/.cache/SatelliteGridding/modis/sinusoidal_2400px_v1}"
CHUNK_DIR="${CHUNK_DIR:-${OUTFILE%.nc}_chunks}"
KEEP_CHUNKS="${KEEP_CHUNKS:-0}"

mkdir -p "$CHUNK_DIR" "$(dirname "$OUTFILE")"

start_epoch=$(date -d "$START_DATE" +%s)
stop_epoch=$(date -d "$STOP_DATE" +%s)
sec_per_day=86400
total_days=$(( (stop_epoch - start_epoch) / sec_per_day + 1 ))
total_windows=$(( (total_days + DDAYS - 1) / DDAYS ))
windows_per_job=$(( (total_windows + JOBS - 1) / JOBS ))
days_per_job=$(( windows_per_job * DDAYS ))

echo "Config:        $CONFIG"
echo "Output:        $OUTFILE"
echo "Window:        $START_DATE -> $STOP_DATE every $DDAYS days"
echo "Resolution:    ${DLAT}deg x ${DLON}deg, global"
echo "Total windows: $total_windows split across $JOBS jobs (~$windows_per_job each)"
echo "Chunk dir:     $CHUNK_DIR"
echo "Geo cache:     $GEO_CACHE"
echo

declare -a CHUNK_FILES=()
declare -a PIDS=()
declare -a JOB_LABELS=()

for ((j=0; j<JOBS; j++)); do
    s_off=$(( j * days_per_job ))
    e_off=$(( (j+1) * days_per_job - 1 ))
    job_start=$(date -I -d "$START_DATE + $s_off days")
    job_stop=$(date  -I -d "$START_DATE + $e_off days")
    [[ "$job_stop" > "$STOP_DATE" ]] && job_stop="$STOP_DATE"
    [[ "$job_start" > "$STOP_DATE" ]] && break

    chunk_idx=$(printf "%02d" "$j")
    chunk_file="$CHUNK_DIR/chunk_${chunk_idx}_${job_start}_${job_stop}.nc"
    log_file="$CHUNK_DIR/chunk_${chunk_idx}.log"
    CHUNK_FILES+=("$chunk_file")
    JOB_LABELS+=("chunk $chunk_idx [$job_start..$job_stop]")

    echo "Launching ${JOB_LABELS[$j]} -> $chunk_file"
    (
        julia --project=. bin/grid.jl center \
            --config "$CONFIG" \
            --latMin -90 --latMax 90 --lonMin -180 --lonMax 180 \
            --dLat "$DLAT" --dLon "$DLON" \
            --startDate "$job_start" --stopDate "$job_stop" \
            --dDays "$DDAYS" \
            --geoProvider modis \
            --geoCache "$GEO_CACHE" \
            --vegIndices \
            --keepGoing \
            -o "$chunk_file" \
            >"$log_file" 2>&1
    ) &
    PIDS+=($!)
done

echo
echo "Waiting for ${#PIDS[@]} jobs..."
fail=0
for ((i=0; i<${#PIDS[@]}; i++)); do
    if wait "${PIDS[$i]}"; then
        echo "  done: ${JOB_LABELS[$i]}"
    else
        echo "  FAILED: ${JOB_LABELS[$i]} (see $CHUNK_DIR/chunk_$(printf "%02d" $i).log)"
        fail=1
    fi
done
[[ "$fail" -eq 0 ]] || { echo "One or more chunks failed."; exit 1; }

echo
echo "Concatenating ${#CHUNK_FILES[@]} chunks -> $OUTFILE"
ncrcat -O -h "${CHUNK_FILES[@]}" "$OUTFILE"
echo "Wrote $OUTFILE"

if [[ "$KEEP_CHUNKS" -eq 0 ]]; then
    rm -f "${CHUNK_FILES[@]}" "$CHUNK_DIR"/chunk_*.log
    rmdir "$CHUNK_DIR" 2>/dev/null || true
fi
