#!/usr/bin/env bash
#
# Regenerate this demo's FDR/TPR figures via the suite-level plotting module
# ../trex_plt_utils.py, using the repo's local .venv. One figure set is
# rendered per base solver (TLARS / TOMP / TAFS); figures are written to
# simulation_results/plots/.
#
# Usage:
#   ./generate_plots.sh                          # all three solver CSVs
#   ./generate_plots.sh <other.csv>              # one explicit CSV
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
    --legend-title 'Dummy distribution' \
    "$@"
fi

for solver in TLARS TOMP TAFS; do
  token="$(echo "$solver" | tr '[:upper:]' '[:lower:]')"
  csv="$here/simulation_results/data/demo_trex_05_dummy_distributions_results_n300_p1000_${token}_random_support.csv"
  "$venv_python" "$here/../trex_plt_utils.py" "$csv" \
    --title 'T-Rex dummy distributions ($n=300$, $p=1000$, random support, 10 runs, '"$solver"' + STANDARD)' \
    --legend-title 'Dummy distribution' \
    "$@"
done
