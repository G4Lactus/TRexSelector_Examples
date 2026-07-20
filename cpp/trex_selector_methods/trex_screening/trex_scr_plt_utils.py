#!/usr/bin/env python3
"""Plotting utilities for the Screen-TRex demonstration suite.

Reads the tidy result CSVs written by the demos (see each demo's
``simulation_results/data/``) and renders publication-style figures into the
sibling ``simulation_results/plots/`` folder. The visual conventions mirror
the T-Rex / DA-TRex / T-Rex+GVS suites: tab10 method colors, hollow per-method
markers, ONE figure-level legend outside the axes, PNG + PDF output.

Two CSV shapes are auto-detected from their columns:

* **Sweep** — ``method,metric,<SNR|rho>,value`` with metrics FDR / TPR /
  Est. FDR (demos 01-03, 06). Rendered as a TPR panel and an FDR panel vs.
  the sweep variable; the FDR panel overlays the procedure's own
  **estimated** FDR as dashed lines so the self-estimate can be compared
  against the realized FDR at a glance.

* **Biobank** — ``method,metric,snr,value`` with an additional ``Usage``
  metric (demos 04-05, "Algorithm 1" adaptive routing). Rendered as three
  panels: TPR, FDR (realized solid / estimated dashed), and Usage % per
  routed method. The FDR panel shades the routing window
  [lower_bound_FDR, upper_bound_FDR] = [0.05, 0.15].

Usage:
    trex_scr_plt_utils.py <results.csv> [--outdir DIR] [--title T]
                          [--formats png pdf] [--stem S]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

# Per-method marker cycle — overlapping curves stay distinguishable.
METHOD_MARKERS = ["o", "s", "^", "D", "v", "P"]
# Routing FDR window of the biobank workflow (BiobankScreenTRexControl).
ROUTING_FDR_WINDOW = (0.05, 0.15)


# ---------------------------------------------------------------------------
# Shared helpers (house conventions)
# ---------------------------------------------------------------------------
def place_figure_legend(fig, handles, labels, group_label: str = "Method",
                        anchor_x: float = 0.82, fontsize: int = 10):
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


def _sweep_column(df: pd.DataFrame) -> str:
    """The sweep-variable column is whatever is not method/metric/value."""
    extras = [c for c in df.columns if c not in ("method", "metric", "value")]
    if len(extras) != 1:
        raise ValueError(f"Expected one sweep column, got: {extras}")
    return extras[0]


def _axis_props(sweep: str) -> tuple[str, bool]:
    """Axis label + log-scale flag for the sweep variable."""
    if sweep.lower() == "snr":
        return "SNR", True
    if sweep.lower() in ("rho", "ρ"):
        return r"$\rho$", False
    return sweep, False


def _plot_series(ax, df, sweep, methods, colors, metric,
                 linestyle="-", markers=True):
    lines = []
    for si, method in enumerate(methods):
        m = df[(df["method"] == method) & (df["metric"] == metric)]
        m = m.sort_values(sweep)
        line, = ax.plot(
            m[sweep], m["value"], linestyle=linestyle,
            marker=METHOD_MARKERS[si % len(METHOD_MARKERS)] if markers else None,
            markersize=5, markerfacecolor="none", markeredgewidth=1.2,
            color=colors[si], label=method)
        lines.append(line)
    return lines


# ---------------------------------------------------------------------------
# Sweep figure (demos 01-03, 06): TPR + FDR/Est.FDR panels
# ---------------------------------------------------------------------------
def plot_sweep(df: pd.DataFrame, title: str) -> plt.Figure:
    sweep = _sweep_column(df)
    xlabel, logx = _axis_props(sweep)
    methods = list(dict.fromkeys(df["method"]))
    colors = plt.get_cmap("tab10").colors

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8))

    # TPR panel
    ax = axes[0]
    tpr_lines = _plot_series(ax, df, sweep, methods, colors, "TPR")
    ax.set_ylim(0, 1.05)
    ax.set_title("TPR (true positive rate)", fontsize=12, fontweight="bold")

    # FDR panel: realized (solid) + estimated (dashed)
    ax = axes[1]
    _plot_series(ax, df, sweep, methods, colors, "FDR")
    _plot_series(ax, df, sweep, methods, colors, "Est. FDR",
                 linestyle="--", markers=False)
    ax.set_title("FDR: realized (solid) vs. estimated (dashed)",
                 fontsize=12, fontweight="bold")

    for ax in axes:
        if logx:
            ax.set_xscale("log")
        ax.grid(alpha=0.3)
        ax.set_xlabel(xlabel, fontsize=11)

    # Legend: method colors, then the linestyle key.
    style_solid = plt.Line2D([], [], color="0.3", linestyle="-")
    style_dash = plt.Line2D([], [], color="0.3", linestyle="--")
    handles = tpr_lines + [style_solid, style_dash]
    labels = methods + ["realized FDR", "estimated FDR"]

    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.98)
    place_figure_legend(fig, handles, labels, "Method",
                        anchor_x=0.815, fontsize=9)
    fig.subplots_adjust(left=0.07, right=0.80, bottom=0.13, top=0.82,
                        wspace=0.22)
    return fig


# ---------------------------------------------------------------------------
# Biobank figure (demos 04-05): TPR + FDR + Usage panels
# ---------------------------------------------------------------------------
def plot_biobank(df: pd.DataFrame, title: str) -> plt.Figure:
    sweep = _sweep_column(df)
    xlabel, logx = _axis_props(sweep)
    methods = list(dict.fromkeys(df["method"]))
    colors = plt.get_cmap("tab10").colors

    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.8))

    ax = axes[0]
    tpr_lines = _plot_series(ax, df, sweep, methods, colors, "TPR")
    ax.set_ylim(0, 1.05)
    ax.set_title("TPR (conditional on routing)", fontsize=12,
                 fontweight="bold")

    ax = axes[1]
    _plot_series(ax, df, sweep, methods, colors, "FDR")
    _plot_series(ax, df, sweep, methods, colors, "Est. FDR",
                 linestyle="--", markers=False)
    lo, hi = ROUTING_FDR_WINDOW
    window = ax.axhspan(lo, hi, color="0.85", alpha=0.6, zorder=0)
    ax.set_title("FDR: realized (solid) vs. estimated (dashed)",
                 fontsize=12, fontweight="bold")

    ax = axes[2]
    usage = df[df["metric"] == "Usage"].copy()
    usage["value"] = usage["value"] * 100.0
    for si, method in enumerate(methods):
        m = usage[usage["method"] == method].sort_values(sweep)
        ax.plot(m[sweep], m["value"],
                marker=METHOD_MARKERS[si % len(METHOD_MARKERS)],
                markersize=5, markerfacecolor="none", markeredgewidth=1.2,
                color=colors[si], label=method)
    ax.set_ylim(0, 105)
    ax.set_title("Usage (% of screenings routed)", fontsize=12,
                 fontweight="bold")

    for ax in axes:
        if logx:
            ax.set_xscale("log")
        ax.grid(alpha=0.3)
        ax.set_xlabel(xlabel, fontsize=11)

    style_solid = plt.Line2D([], [], color="0.3", linestyle="-")
    style_dash = plt.Line2D([], [], color="0.3", linestyle="--")
    handles = tpr_lines + [style_solid, style_dash, window]
    labels = methods + ["realized FDR", "estimated FDR",
                        f"routing window [{lo:g}, {hi:g}]"]

    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.99)
    place_figure_legend(fig, handles, labels, "Method",
                        anchor_x=0.845, fontsize=9)
    fig.subplots_adjust(left=0.05, right=0.83, bottom=0.13, top=0.83,
                        wspace=0.24)
    return fig


# ---------------------------------------------------------------------------
# Dispatch + CLI
# ---------------------------------------------------------------------------
def detect_kind(df: pd.DataFrame) -> str:
    if not {"method", "metric", "value"} <= set(df.columns):
        raise ValueError(f"Unrecognized CSV schema: {sorted(df.columns)}")
    if "Usage" in set(df["metric"]):
        return "biobank"
    return "sweep"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Render Screen-TRex result figures.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", type=Path, help="Result CSV (schema auto-detected).")
    p.add_argument("--outdir", type=Path, default=None,
                   help="Output dir (default: sibling simulation_results/plots).")
    p.add_argument("--title", default=None, help="Figure title.")
    p.add_argument("--formats", nargs="+", default=["png", "pdf"])
    p.add_argument("--stem", default=None,
                   help="Output filename stem (default: CSV stem).")
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

    kind = detect_kind(df)
    outdir = args.outdir or default_plots_dir(args.csv)
    stem = args.stem or args.csv.stem
    title = args.title or args.csv.stem

    fig = plot_biobank(df, title) if kind == "biobank" else plot_sweep(df, title)
    save_figure(fig, outdir, stem, args.formats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
