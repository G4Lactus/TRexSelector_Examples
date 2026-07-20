#!/usr/bin/env python3
"""Plotting utilities for the T-Rex SPCA demonstration suite.

Reads the tidy result CSV written by the demo (see
``demo_trex_spca_01_mc_sim/simulation_results/data/``) and renders a
publication-style figure into the sibling ``simulation_results/plots/`` folder.
The visual conventions mirror the other suites: tab10 method colours, hollow
per-method markers, ONE figure-level legend outside the axes, PNG + PDF output.

CSV shape (``method,metric,snr_db,value``): the sweep axis is the
signal-to-noise ratio **in decibels**, and it is plotted as such — the grid is
linear in dB (-10 … 10) and is never converted back to a linear SNR. Three
metrics are reported per method:

* **FDR** — false discovery rate of the estimated loading support (PC1),
  against the ``tfdr`` target drawn as a dashed reference line.
* **TPR** — true positive rate of that support.
* **PEV** — proportion of explained variance, the quantity the sparsity is
  bought against. PCA baselines lead here by construction, so the panel shows
  what the FDR control costs.

Usage:
    trex_spca_plt_utils.py <results.csv> [--outdir DIR] [--title T]
                           [--tfdr F] [--formats png pdf] [--stem S]
                           [--xscale {linear,sqrt}]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import (  # noqa: E402
    AutoMinorLocator, MaxNLocator, ScalarFormatter,
)

DEFAULT_TFDR = 0.1
# Per-method marker cycle — overlapping curves stay distinguishable.
METHOD_MARKERS = ["o", "s", "^", "D", "v", "P"]
# Baselines first, then the T-Rex SPCA variants, so the legend reads top-down.
METHOD_ORDER = ["OrdPCA", "OraclePCA",
                "TRexSPCA-EN-Act", "TRexSPCA-EN-Thr",
                "TRexSPCA-ENaug-Act", "TRexSPCA-ENaug-Thr"]
METRICS = [("TPR", "TPR (true positive rate)", (0, 1.05)),
           ("FDR", "FDR (false discovery rate)", None),
           ("PEV", "PEV (proportion of explained variance)", None)]


# ---------------------------------------------------------------------------
# Shared helpers (house conventions)
# ---------------------------------------------------------------------------
def place_figure_legend(fig, handles, labels, group_label: str = "Method",
                        anchor_x: float = 0.845, fontsize: int = 9):
    """One figure-level legend outside the right panel, opaque white box."""
    header = plt.Line2D([], [], linestyle="None", alpha=0.0)
    handles = [header] + handles
    labels = [f"{group_label}:"] + labels
    legend = fig.legend(
        handles, labels,
        loc="center left", bbox_to_anchor=(anchor_x, 0.5),
        frameon=True, framealpha=1.0, facecolor="white", edgecolor="0.7",
        fontsize=fontsize, ncol=1,
    )
    legend.get_frame().set_linewidth(0.8)
    return legend


def apply_db_xaxis(ax, values: list[float], xscale: str = "linear") -> None:
    """Configure the decibel sweep axis.

    The grid is already logarithmic in nature — decibels *are* a log scale — so
    it is plotted linearly in dB and never re-mapped. Major ticks come from the
    locators rather than from the swept values, so labels cannot collide, and
    minor ticks carry a faint grid. ``sqrt`` is offered for parity with the
    other suites but is rarely useful on an axis that already spans evenly.
    """
    lo, hi = min(values), max(values)
    span = hi - lo
    pad = span * 0.04 if span > 0 else 0.5
    if xscale == "sqrt":
        shift = -lo + 1.0 if lo <= 0 else 0.0    # keep the transform defined
        ax.set_xscale("function",
                      functions=(lambda x: np.sqrt(np.maximum(x + shift, 0.0)),
                                 lambda x: np.square(x) - shift))
    ax.set_xlim(lo - pad, hi + pad)
    ax.xaxis.set_major_locator(MaxNLocator(nbins="auto", steps=[1, 2, 2.5, 5, 10]))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.grid(which="minor", axis="x", alpha=0.25, linewidth=0.5)


def save_figure(fig: plt.Figure, outdir: Path, stem: str,
                formats: list[str]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        path = outdir / f"{stem}.{ext}"
        fig.savefig(path, dpi=300 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
        print(f"Wrote: {path}")
    plt.close(fig)


def default_plots_dir(csv_path: Path) -> Path:
    """`.../simulation_results/data/x.csv` -> `.../simulation_results/plots`."""
    if csv_path.parent.name == "data":
        return csv_path.parent.parent / "plots"
    return csv_path.parent / "plots"


# ---------------------------------------------------------------------------
# Figure: TPR / FDR / PEV vs. SNR (dB)
# ---------------------------------------------------------------------------
def plot_spca_sweep(df: pd.DataFrame, title: str, tfdr: float = DEFAULT_TFDR,
                    xscale: str = "linear") -> plt.Figure:
    present = list(dict.fromkeys(df["method"]))
    methods = [m for m in METHOD_ORDER if m in present]
    methods += [m for m in present if m not in methods]   # keep any newcomers
    xs = sorted(df["snr_db"].unique())
    colors = plt.get_cmap("tab10").colors

    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.8))
    legend_handles: list = []
    legend_labels: list[str] = []

    for ai, (metric, mlabel, ylim) in enumerate(METRICS):
        ax = axes[ai]
        for mi, method in enumerate(methods):
            sub = df[(df["method"] == method) & (df["metric"] == metric)]
            sub = sub.sort_values("snr_db")
            line, = ax.plot(
                sub["snr_db"], sub["value"],
                marker=METHOD_MARKERS[mi % len(METHOD_MARKERS)],
                markersize=5, markerfacecolor="none", markeredgewidth=1.2,
                color=colors[mi], label=method)
            if ai == 0:
                legend_handles.append(line)
                legend_labels.append(method)
        if ylim:
            ax.set_ylim(*ylim)
        if metric == "FDR" and tfdr > 0:
            ref = ax.axhline(tfdr, color="red", linestyle="--", linewidth=1.2)
            legend_handles.append(ref)
            legend_labels.append(f"tFDR = {tfdr:g}")
        ax.set_title(mlabel, fontsize=12, fontweight="bold")
        ax.grid(alpha=0.3)
        apply_db_xaxis(ax, xs, xscale)
        ax.set_xlabel("SNR (dB)", fontsize=11)

    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.99)
    place_figure_legend(fig, legend_handles, legend_labels, "Method")
    fig.subplots_adjust(left=0.05, right=0.83, bottom=0.13, top=0.83,
                        wspace=0.24)
    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Render T-Rex SPCA result figures.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", type=Path, help="Result CSV (method,metric,snr_db,value).")
    p.add_argument("--outdir", type=Path, default=None,
                   help="Output dir (default: sibling simulation_results/plots).")
    p.add_argument("--title", default=None, help="Figure title.")
    p.add_argument("--tfdr", type=float, default=DEFAULT_TFDR)
    p.add_argument("--formats", nargs="+", default=["png", "pdf"])
    p.add_argument("--stem", default=None,
                   help="Output filename stem (default: CSV stem).")
    p.add_argument("--xscale", choices=["linear", "sqrt"], default="linear",
                   help="Sweep-axis scale; the axis is in dB either way.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {args.csv}")
    df = pd.read_csv(args.csv, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    for c in df.columns:
        if not pd.api.types.is_numeric_dtype(df[c]):
            df[c] = df[c].astype(str).str.strip()

    missing = {"method", "metric", "snr_db", "value"} - set(df.columns)
    if missing:
        raise ValueError(f"CSV is missing columns: {sorted(missing)}")

    fig = plot_spca_sweep(df, args.title or args.csv.stem, args.tfdr,
                          args.xscale)
    save_figure(fig, args.outdir or default_plots_dir(args.csv),
                args.stem or args.csv.stem, args.formats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
