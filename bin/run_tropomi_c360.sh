#!/usr/bin/env bash
#
# Grid TROPOMI Level-2 SIF onto a GEOS/GCHP cubed-sphere grid (default C360) as a
# 16-day product. Same fan-out shape as run_tropomi_l2_16day.sh: independent
# Julia processes over time chunks, then `ncrcat` along the time dimension.
#
# A cubed sphere is global, so there are no LAT_MIN/LON_MIN-style bounds — only
# the resolution Nc. Gridding uses the sequential CPU path (the cube does not yet
# have a GPU/KA path); parallelism comes from the per-chunk processes.
#
# Run from the SatelliteGridding repo root:
#   CONFIG=examples/tropomi_sif.toml NC=360 bash bin/run_tropomi_c360.sh
#
# Tunables (env vars):
#   CONFIG=examples/tropomi_sif.toml    L2 config (must define basic.lat_bnd/lon_bnd)
#   NC=360                              cubed-sphere cells per panel edge (C<NC>)
#   CS_CONVENTION=gmao                  gmao (GEOS/GCHP, -10° offset) or equiangular
#   JOBS=8                              concurrent Julia processes
#   OUTFILE=...                         final concatenated netCDF
#   CHUNK_DIR=...                       per-chunk intermediates dir
#   START_DATE=2018-05-01
#   STOP_DATE=2025-12-31
#   DDAYS=16
#   N_OVERSAMPLE=0                      0 means the cubed-sphere default (4)
#   KEEP_CHUNKS=0
#   PROGRESS_INTERVAL=30

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

CONFIG="${CONFIG:-examples/tropomi_sif.toml}"
NC="${NC:-360}"
CS_CONVENTION="${CS_CONVENTION:-gmao}"
PRODUCT_TAG="${PRODUCT_TAG:-$(basename "${CONFIG%.toml}")}"
OUTFILE="${OUTFILE:-${PRODUCT_TAG}_16d_c${NC}.nc}"
START_DATE="${START_DATE:-2018-05-01}"
STOP_DATE="${STOP_DATE:-2025-12-31}"
DDAYS="${DDAYS:-16}"
JOBS="${JOBS:-8}"
N_OVERSAMPLE="${N_OVERSAMPLE:-0}"
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
echo "Grid:          cubed sphere C${NC} ($CS_CONVENTION convention)"
echo "Window:        $START_DATE -> $STOP_DATE every $DDAYS days"
echo "Oversample:    n_oversample=$N_OVERSAMPLE (0=cubed-sphere default)"
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
            --gridType cs --Nc "$NC" --csConvention "$CS_CONVENTION" \
            --startDate "$job_start" --stopDate "$job_stop" \
            --dDays "$DDAYS" \
            --nOversample "$N_OVERSAMPLE" \
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
set +e
set +o pipefail
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
set -e
set -o pipefail

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
