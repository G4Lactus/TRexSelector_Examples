#!/usr/bin/env bash
#
# Regenerate this demo's FDR/TPR figures via the plotting module
#   ../trex_plt_utils.py
# using the repo's local .venv. Figures are written to
# simulation_results/plots/.
#
# Usage:
#   ./generate_plots.sh                          # this demo's default CSV
#   ./generate_plots.sh <other.csv>              # explicit CSV
#   ./generate_plots.sh --formats png            # extra args pass through to
#   ./generate_plots.sh --tfdr 0.05              # trex_plt_utils.py
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Find the repo-root .venv by walking up the directory tree.
dir="$here"
venv_python=""
while [[ "$dir" != "/" ]]; do
  if [[ -x "$dir/.venv/bin/python" ]]; then
    venv_python="$dir/.venv/bin/python"
    break
  fi
  dir="$(dirname "$dir")"
done

if [[ -z "$venv_python" ]]; then
  echo "error: no .venv/bin/python found above $here" >&2
  echo "       create it and 'pip install matplotlib pandas plotly', then re-run." >&2
  exit 1
fi

# First non-flag argument: plot just that CSV instead of the three defaults.
if [[ $# -gt 0 && "$1" != -* ]]; then
  csv="$1"
  shift
  exec "$venv_python" "$here/../trex_plt_utils.py" "$csv" \
    --legend-title 'L-strategy' \
    "$@"
fi

# TAFS is the only base solver here run with a non-zero AFS correlation
# parameter (rho_afs = 0.3 in this demo's solver table; TLARS/TOMP take 0), so
# the figures name it -- otherwise the setting is lost when a plot is read apart
# from this demo's README. The loop variable stays the raw name: it also builds
# the CSV filename token.
solver_label() {
  case "$1" in
    TAFS) echo 'TAFS (rho = 0.3)' ;;
    *)    echo "$1" ;;
  esac
}

# One per-solver figure set (overview + grouped + interactive html) per base
# solver, via the suite-level plotter.
for solver in TLARS TOMP TAFS; do
  token="$(echo "$solver" | tr '[:upper:]' '[:lower:]')"
  label="$(solver_label "$solver")"
  csv="$here/simulation_results/data/demo_trex_04_lloop_strategies_results_n300_p1000_${token}_random_support.csv"
  "$venv_python" "$here/../trex_plt_utils.py" "$csv" \
    --title 'T-Rex L-loop strategies ($n=300$, $p=1000$, random support, '"$label"')' \
    --legend-title 'L-strategy' \
    "$@"
done

# Cross-solver comparison grid (2x3: FDR/TPR rows x TLARS/TOMP/TAFS cols, all
# eight L-loop strategies per panel on a shared scale) -- the suite plotter's
# 'grid' mode. This is the figure for the demo's core question: whether the
# SKIPL FDR overshoot is specific to the LARS path (TLARS) or persists under the
# greedy solvers (TOMP / TAFS).
data="$here/simulation_results/data"
"$venv_python" "$here/../trex_plt_utils.py" grid \
  --labels TLARS TOMP "$(solver_label TAFS)" \
  --csvs \
    "$data/demo_trex_04_lloop_strategies_results_n300_p1000_tlars_random_support.csv" \
    "$data/demo_trex_04_lloop_strategies_results_n300_p1000_tomp_random_support.csv" \
    "$data/demo_trex_04_lloop_strategies_results_n300_p1000_tafs_random_support.csv" \
  --title 'T-Rex L-loop strategies across base solvers ($n=300$, $p=1000$, random support)' \
  --legend-title 'L-strategy' \
  --stem demo_trex_04_lloop_strategies_grid \
  --outdir "$here/simulation_results/plots" \
  "$@"
