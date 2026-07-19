#!/bin/bash
# ==============================================================================
# run_gvs_demos_overnight.sh
# ==============================================================================
# Run the trex_gvs demos sequentially as an unattended (overnight) job.
#
# Default order: 02 03 04 05 06 07 08  (demo 01 is excluded because it is
# typically launched separately; pass demo numbers to override, e.g.
#   ./run_gvs_demos_overnight.sh 01 02 03
# runs exactly those).
#
# Behaviour:
#   - rebuilds the selected demo targets first (cmake build/release)
#   - runs each demo binary in order; a failure does not stop the chain
#   - per-demo console logs in logs/overnight_<timestamp>/demo_XX.log
#     (result TXT/CSV files go to each demo's own simulation_results/data/,
#      written by the demos themselves)
#   - keeps the Mac awake while running (caffeinate -i, tied to this script)
#   - prints a summary table with per-demo exit codes and wall-clock times
#
# Launch for an overnight run:
#   cd cpp/trex_selector_methods/trex_gvs
#   nohup ./run_gvs_demos_overnight.sh > /dev/null 2>&1 &
# then follow progress with:  tail -f logs/overnight_*/demo_02*.log
# ==============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPP_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BIN_DIR="${CPP_ROOT}/build/release/bin/trex_selector_methods/trex_gvs"

# Demo number -> full name
demo_name() {
    case "$1" in
        01) echo "demo_trex_gvs_01_mc_sim_hastie_en_blocks" ;;
        02) echo "demo_trex_gvs_02_mc_sim_scattered_blocks" ;;
        03) echo "demo_trex_gvs_03_mc_sim_mixed_blocks" ;;
        04) echo "demo_trex_gvs_04_mc_sim_neg_traps" ;;
        05) echo "demo_trex_gvs_05_mc_sim_t3_blocks" ;;
        06) echo "demo_trex_gvs_06_mc_sim_ar1_blocks" ;;
        07) echo "demo_trex_gvs_07_mc_sim_arma_blocks" ;;
        08) echo "demo_trex_gvs_08_mc_sim_block_bench" ;;
        *)  echo "" ;;
    esac
}

DEMOS=("$@")
if [ ${#DEMOS[@]} -eq 0 ]; then
    DEMOS=(02 03 04 05 06 07 08)
fi

TS="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${SCRIPT_DIR}/logs/overnight_${TS}"
mkdir -p "${LOG_DIR}"

echo "trex_gvs overnight run — $(date)"
echo "Demos: ${DEMOS[*]}"
echo "Logs:  ${LOG_DIR}"
echo

# Keep the machine awake for the lifetime of this script (macOS).
if command -v caffeinate >/dev/null 2>&1; then
    caffeinate -i -w $$ &
fi

# ------------------------------------------------------------------------------
# Rebuild the selected targets so the binaries match the current sources.
# ------------------------------------------------------------------------------
TARGETS=()
for d in "${DEMOS[@]}"; do
    name="$(demo_name "$d")"
    if [ -z "$name" ]; then
        echo "Unknown demo number: $d — skipping."
        continue
    fi
    TARGETS+=("$name")
done

echo "Building: ${TARGETS[*]}"
if ! cmake --build "${CPP_ROOT}/build/release" --target "${TARGETS[@]}" -j 8 \
        > "${LOG_DIR}/build.log" 2>&1; then
    echo "BUILD FAILED — see ${LOG_DIR}/build.log. Aborting."
    exit 1
fi
echo "Build OK."
echo

# ------------------------------------------------------------------------------
# Run demos sequentially.
# ------------------------------------------------------------------------------
declare -a SUMMARY

for d in "${DEMOS[@]}"; do
    name="$(demo_name "$d")"
    [ -z "$name" ] && continue
    bin="${BIN_DIR}/${name}/${name}"
    log="${LOG_DIR}/demo_${d}.log"

    if [ ! -x "$bin" ]; then
        echo "[$(date +%H:%M:%S)] demo ${d}: binary not found (${bin}) — SKIPPED"
        SUMMARY+=("demo ${d}  MISSING")
        continue
    fi

    echo "[$(date +%H:%M:%S)] demo ${d} (${name}) starting ..."
    t0=$(date +%s)
    "$bin" > "$log" 2>&1
    rc=$?
    t1=$(date +%s)
    dt=$((t1 - t0))
    hms=$(printf '%dh %02dm %02ds' $((dt/3600)) $(((dt%3600)/60)) $((dt%60)))

    if [ $rc -eq 0 ]; then
        echo "[$(date +%H:%M:%S)] demo ${d} finished OK in ${hms}"
        SUMMARY+=("demo ${d}  OK      ${hms}")
    else
        echo "[$(date +%H:%M:%S)] demo ${d} FAILED (exit ${rc}) after ${hms} — see ${log}"
        SUMMARY+=("demo ${d}  FAILED(${rc})  ${hms}")
    fi
done

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo
echo "=============================================================="
echo "  trex_gvs overnight run summary — $(date)"
echo "=============================================================="
for line in "${SUMMARY[@]}"; do
    echo "  $line"
done
echo "  Logs: ${LOG_DIR}"
echo "=============================================================="
