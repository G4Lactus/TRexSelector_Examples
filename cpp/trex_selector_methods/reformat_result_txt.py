#!/usr/bin/env python3
"""Rewrite saved MC result tables in the current console layout.

The demos print (and save) their Monte Carlo summary tables from the
``save_and_print_*`` helpers in each suite's simulation-utils header. When that
layout changes, the ``.txt`` files already saved under each demo's
``simulation_results/data/`` go stale. Rather than re-running hours of Monte
Carlo, this script regenerates them from the sibling ``.csv``, which carries the
same numbers in tidy long format.

It mirrors the C++ renderer: the series name on its own line, metric rows
indented beneath it, and the sweep axis as a single header row. Any framed
metadata banner at the top of the old ``.txt`` (the DA-TRex "=== ... MC=200 ..."
block) is preserved verbatim.

Only the columnar ``<series>/<metric>/<sweep>`` tables are rewritten. Files
whose CSV holds a two-dimensional matrix (an extra label column, e.g.
``row_label``/``col_label``, ``support_type``, ``gap``) use a different
renderer and are skipped.

Caveat: the CSV stores 6 decimals while the table shows 4, so a value sitting
exactly on a rounding tie can differ by one unit in the last displayed digit
from a fresh C++ run. Re-run the demo if bit-exact tables matter.

Usage:
    ./reformat_result_txt.py                 # every suite below this folder
    ./reformat_result_txt.py trex trex_da    # only these suites
    ./reformat_result_txt.py <file.csv> ...  # only these files
"""
from __future__ import annotations

import csv
import re
import sys
from collections import OrderedDict
from pathlib import Path

INDENT_W = 2
METRIC_W = 23
COL_W = 10

# CSV sweep column -> sweep label as printed by the C++ renderer.
SWEEP_LABELS = {"snr": "SNR", "snr_db": "SNR(dB)", "rho": "Rho"}

# CSV metric tag -> label as printed in the table.
METRIC_LABELS = {"AvgL": "Avg L", "AvgT": "Avg T",
                 "Avg_L": "Avg L", "Avg_T": "Avg T"}

# A tidy columnar CSV has exactly these four columns (in any order):
#   <series>, metric, <sweep>, value
SERIES_KEYS = ("solver", "method")


def classify(path: Path):
    """Return (series_col, sweep_col) for a columnar CSV, or None if unsupported."""
    with open(path, newline="") as fh:
        header = next(csv.reader(fh), [])
    if len(header) != 4 or "metric" not in header or "value" not in header:
        return None
    series = next((c for c in header if c in SERIES_KEYS), None)
    if series is None:
        return None
    sweep = next(c for c in header
                 if c not in (series, "metric", "value"))
    return series, sweep


def read_table(path: Path, series_col: str, sweep_col: str):
    rows = list(csv.DictReader(open(path)))
    xs, series, metrics = [], OrderedDict(), OrderedDict()
    values = {}
    for r in rows:
        x = float(r[sweep_col])
        if x not in xs:
            xs.append(x)
        series.setdefault(r[series_col], None)
        metrics.setdefault(r["metric"], None)
        values[(r[series_col], r["metric"], x)] = float(r["value"])
    return sorted(xs), list(series), list(metrics), values


def banner(txt_path: Path) -> str:
    """Everything above the old table header (framed metadata), verbatim."""
    if not txt_path.exists():
        return ""
    text = txt_path.read_text()
    m = re.search(r"^\s*(Solver|Method)\s+Metric\b", text, re.M)
    if m:                       # old columnar layout
        return text[:m.start()].rstrip("\n") + "\n"
    m = re.search(r"^\s{2}\S.*$", text, re.M)   # already-new layout: keep header block
    parts = text.split("\n\n")
    return parts[0].rstrip("\n") + "\n" if parts else ""


def render(head: str, sweep_label: str, xs, series, metrics, values) -> str:
    out = [head.rstrip("\n"), ""]
    hdr = " " * INDENT_W + sweep_label.ljust(METRIC_W)
    hdr += "".join(f"{x:>{COL_W}.2f}" for x in xs)
    out.append(hdr)
    out.append("-" * (INDENT_W + METRIC_W + COL_W * len(xs)))
    for s in series:
        out.append("")
        out.append(s)
        for m in metrics:
            if (s, m, xs[0]) not in values:
                continue
            label = METRIC_LABELS.get(m, m)
            row = " " * INDENT_W + label.ljust(METRIC_W)
            row += "".join(f"{values[(s, m, x)]:>{COL_W}.4f}" for x in xs)
            out.append(row)
    out.append("")
    return "\n".join(out) + "\n"


def main(argv: list[str]) -> int:
    here = Path(__file__).resolve().parent
    args = argv[1:]
    if args and all(a.endswith(".csv") for a in args):
        csvs = [Path(a) for a in args]
    else:
        suites = args or ["trex", "trex_da", "trex_gvs", "trex_screening", "trex_spca"]
        csvs = sorted(c for s in suites
                      for c in (here / s).glob("demo_*/simulation_results/data/*.csv"))

    done = skipped = 0
    for c in csvs:
        kind = classify(c)
        if kind is None:
            skipped += 1
            continue
        series_col, sweep_col = kind
        xs, series, metrics, values = read_table(c, series_col, sweep_col)
        txt = c.with_suffix(".txt")
        if not txt.exists():
            # No sibling table: writing one would drop the framed metadata
            # banner that only the saved .txt carries.
            print(f"Skipped (no .txt): {c.name}")
            skipped += 1
            continue
        txt.write_text(render(banner(txt),
                              SWEEP_LABELS.get(sweep_col, sweep_col),
                              xs, series, metrics, values))
        try:
            shown = txt.relative_to(here)
        except ValueError:
            shown = txt
        print(f"Rewrote: {shown}")
        done += 1
    print(f"\n{done} table(s) rewritten, {skipped} skipped (matrix or unsupported schema)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
