#!/usr/bin/env bash
#
# Grid TROPOMI Level-2 data into a 16-day product, defaulting to 1° global.
# Override LAT_MIN/LAT_MAX/LON_MIN/LON_MAX for a subregion. Footprints come
# from `lat_bnd`/`lon_bnd` in the L2 file, so this uses the `l2` subcommand
# of `grid.jl` (subpixel oversampling, not center-coord binning).
#
# Same overall shape as run_modis_16day_global.sh: fans out across time
# chunks as independent Julia processes, then `ncrcat`s the per-chunk netCDFs
# into one output along the time dimension.
#
# Run from the SatelliteGridding repo root:
#   CONFIG=examples/tropomi_sif.toml bash bin/run_tropomi_l2_16day.sh
#
# Tunables (env vars):
#   CONFIG=examples/tropomi_sif.toml    L2 config (must define basic.lat_bnd/lon_bnd)
#   JOBS=8                              concurrent Julia processes
#   OUTFILE=...                         final concatenated netCDF
#   CHUNK_DIR=...                       per-chunk intermediates dir
#   START_DATE=2018-05-01
#   STOP_DATE=2025-12-31
#   DDAYS=16
#   DLAT=1.0   DLON=1.0
#   LAT_MIN=-90 LAT_MAX=90 LON_MIN=-180 LON_MAX=180
#   N_OVERSAMPLE=0                      0 means auto-compute from footprint/grid ratio
#   FOOTPRINT=quad                      or "circle" for circular L2 products
#   KEEP_CHUNKS=1                       keep per-chunk files after concat
#   PROGRESS_INTERVAL=30                seconds between live progress printouts

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

CONFIG="${CONFIG:-examples/tropomi_sif.toml}"
PRODUCT_TAG="${PRODUCT_TAG:-$(basename "${CONFIG%.toml}")}"
OUTFILE="${OUTFILE:-${PRODUCT_TAG}_16d_1deg.nc}"
START_DATE="${START_DATE:-2018-05-01}"
STOP_DATE="${STOP_DATE:-2025-12-31}"
DDAYS="${DDAYS:-16}"
DLAT="${DLAT:-1.0}"
DLON="${DLON:-1.0}"
LAT_MIN="${LAT_MIN:--90}"
LAT_MAX="${LAT_MAX:-90}"
LON_MIN="${LON_MIN:--180}"
LON_MAX="${LON_MAX:-180}"
JOBS="${JOBS:-8}"
N_OVERSAMPLE="${N_OVERSAMPLE:-0}"
FOOTPRINT="${FOOTPRINT:-quad}"
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
echo "Resolution:    ${DLAT}deg x ${DLON}deg, lat[$LAT_MIN, $LAT_MAX] lon[$LON_MIN, $LON_MAX]"
echo "Footprint:     $FOOTPRINT (n_oversample=$N_OVERSAMPLE; 0=auto)"
echo "Total windows: $total_windows split across $JOBS jobs (~$windows_per_job each)"
echo "Chunk dir:     $CHUNK_DIR"
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
        julia --project=. bin/grid.jl l2 \
            --config "$CONFIG" \
            --latMin "$LAT_MIN" --latMax "$LAT_MAX" \
            --lonMin "$LON_MIN" --lonMax "$LON_MAX" \
            --dLat "$DLAT" --dLon "$DLON" \
            --startDate "$job_start" --stopDate "$job_stop" \
            --dDays "$DDAYS" \
            --nOversample "$N_OVERSAMPLE" \
            --footprint "$FOOTPRINT" \
            --keepGoing \
            -o "$chunk_file" \
            >"$log_file" 2>&1
    ) &
    PIDS+=($!)
done

echo
echo "Waiting for ${#PIDS[@]} jobs..."
echo "  (live tail: tail -F $CHUNK_DIR/chunk_*.log)"

PROGRESS_INTERVAL="${PROGRESS_INTERVAL:-30}"
while :; do
    any_running=0
    line=""
    for ((i=0; i<${#PIDS[@]}; i++)); do
        if kill -0 "${PIDS[$i]}" 2>/dev/null; then
            any_running=1
            log="$CHUNK_DIR/chunk_$(printf "%02d" $i).log"
            cleaned=$(tail -c 2000 "$log" 2>/dev/null \
                | tr -d '\r' | sed 's/\x1b\[[A-Za-z]//g; s/\[[A-Za-z]//g')
            pct=$(echo "$cleaned" | grep -oE "Files: +[0-9]+%" | tail -1)
            eta=$(echo "$cleaned" | grep -oE "ETA: +[0-9:]+" | tail -1)
            line+="    chunk $(printf "%02d" $i): ${pct:-starting...}  ${eta}"$'\n'
        fi
    done
    [[ "$any_running" -eq 0 ]] && break
    printf "\n[%s]\n%s" "$(date '+%H:%M:%S')" "$line"
    sleep "$PROGRESS_INTERVAL"
done

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
