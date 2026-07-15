#!/usr/bin/env bash
#
# Regenerate this demo's scalability dashboard via the suite plotter's
# 'scalability' mode:
#   ../trex_plt_utils.py scalability <csv> ...
# using the repo's local .venv. Figures are written to
# simulation_results/plots/.
#
# Demo 08's CSV uses a wider schema (scenario,solver,n,p,num_mc,snr,metric,
# value) than the rest of the suite; the shared module reads it via
# read_scalability() and renders it with plot_scalability_dashboard().
#
# Usage:
#   ./generate_plots.sh                          # this demo's default CSV
#   ./generate_plots.sh <other.csv>              # explicit CSV
#   ./generate_plots.sh --formats png            # extra args pass through to
#   ./generate_plots.sh --tfdr 0.05              # the scalability plotter
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
  echo "       create it and 'pip install matplotlib pandas', then re-run." >&2
  exit 1
fi

# First non-flag argument overrides this demo's default (wide-schema) CSV.
csv="$here/simulation_results/data/demo_trex_08_scalability_results.csv"
if [[ $# -gt 0 && "$1" != -* ]]; then
  csv="$1"
  shift
fi

# The scalability dashboard lives in the suite plotter's 'scalability' mode.
exec "$venv_python" "$here/../trex_plt_utils.py" scalability "$csv" \
  --title 'T-Rex scalability on the memory-mapped pipeline (SNR grid $\{0.5, 1.0, 2.0\}$)' \
  --stem demo_trex_08_scalability \
  --outdir "$here/simulation_results/plots" \
  "$@"
