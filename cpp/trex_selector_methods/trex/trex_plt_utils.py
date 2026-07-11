#!/usr/bin/env python3
"""
trex_plt_utils.py — shared plotting utilities for the classical T-Rex demo suite.

Suite-level companion of trex_sim_utils.hpp: one module hosting the plot
patterns for every demo in this folder whose results use the tidy CSV schema

    solver,metric,snr,value

(the "solver" column holds whatever the demo sweeps: base solvers in demos
02/03, L-loop strategies in demo 04, dummy distributions in demo 05, ...).
Row colours are inferred from the data, so no per-demo configuration is
required; each demo's generate_plots.sh supplies its CSV path and titles.

Plot patterns hosted here:
    plot_fdr_tpr_vs_snr()          — side-by-side FDR | TPR panels (static).
    plot_fdr_tpr_vs_snr_grouped()  — de-cluttered 2x2 grid, rows split in half.
    write_fdr_tpr_vs_snr_html()    — interactive Plotly HTML (hover, toggles).

Usage (typically via a demo's generate_plots.sh):
    python trex_plt_utils.py <path/to/results.csv>
    python trex_plt_utils.py <csv> --title "..." --legend-title "L-strategy" \
        --stem my_plot --formats png pdf

Output location: figures default to the simulation_results/plots/ sibling when
the CSV lives in simulation_results/data/ (the suite convention); otherwise
they are written next to the CSV. Override with --outdir.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless / batch-safe: no display required

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import ScalarFormatter

# ---------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------
STEM_SUFFIX = "_fdr_tpr_vs_snr"
DEFAULT_TFDR = 0.10          # target FDR control level (tFDR)
DEFAULT_LEGEND_TITLE = "Solver"
Y_HEADROOM = 1.25            # budget above the largest FDR value / the tFDR line
MAX_XTICKS = 6               # thin x-axis ticks when the SNR grid is dense

METRICS = ("FDR", "TPR")
METRIC_LABELS = {
    "FDR": "False discovery rate",
    "TPR": "True positive rate",
}

# Black is reserved for the tFDR reference line, so no data row uses it.
TFDR_COLOR = "black"

# Markers cycle independently of colour so many rows stay distinguishable.
MARKERS = [
    "o", "s", "^", "D", "v", "P", "X",
    "<", ">", "h", "H", "*", "8", "p",
]


def read_mc_results(path: Path) -> pd.DataFrame:
    """Load and validate a tidy Monte Carlo-results CSV (solver,metric,snr,value)."""
    if not path.exists():
        raise FileNotFoundError(
            f"Input file not found:\n  {path}\n"
            "Pass the results CSV path as the first argument."
        )

    df = pd.read_csv(path)
    df.columns = df.columns.str.strip().str.lower()

    required = {"solver", "metric", "snr", "value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"CSV is missing required columns: {sorted(missing)}\n"
            f"Found columns: {df.columns.tolist()}"
        )

    df = df[["solver", "metric", "snr", "value"]].copy()
    df["solver"] = df["solver"].astype(str).str.strip()
    df["metric"] = df["metric"].astype(str).str.strip().str.upper()
    df["snr"] = pd.to_numeric(df["snr"], errors="raise")
    df["value"] = pd.to_numeric(df["value"], errors="raise")

    present = set(df["metric"])
    missing_metrics = set(METRICS) - present
    if missing_metrics:
        raise ValueError(
            f"CSV is missing required metric(s): {sorted(missing_metrics)}\n"
            f"Found metrics: {sorted(present)}"
        )

    # Some demos also record extra metrics (e.g. AvgL, AvgT); keep only the
    # metrics this module plots and note the rest.
    extra = present - set(METRICS)
    if extra:
        print(f"note: ignoring non-plotted metric(s): {sorted(extra)}")
    df = df.loc[df["metric"].isin(METRICS)]

    return df.sort_values(["metric", "solver", "snr"])


def default_plots_dir(csv_path: Path) -> Path:
    """Suite convention: CSVs live in simulation_results/data/, figures in the
    sibling simulation_results/plots/. Falls back to the CSV's own directory
    for CSVs outside that layout."""
    if csv_path.parent.name == "data":
        return csv_path.parent.parent / "plots"
    return csv_path.parent


def row_colors(rows: list[str]) -> dict[str, tuple]:
    """Assign a stable, clearly distinct colour to each data row.

    Concatenates two *qualitative* colormaps (``tab10`` + ``Dark2`` = 18 hand-
    picked, high-contrast colours) rather than sampling a smooth gradient, so
    adjacent curves are easy to tell apart. Colours are keyed by the row's
    position in the (file-order) list, so the mapping is deterministic and works
    for any demo's row set without a hardcoded table.
    """
    palette = list(plt.get_cmap("tab10").colors) + list(
        plt.get_cmap("Dark2").colors
    )
    return {
        row: palette[index % len(palette)]
        for index, row in enumerate(rows)
    }


def rgba_to_css(color: tuple) -> str:
    """Convert a matplotlib RGB(A) tuple in [0, 1] to a CSS ``rgb(...)`` string."""
    r, g, b = (int(round(channel * 255)) for channel in color[:3])
    return f"rgb({r},{g},{b})"


def split_rows(rows: list[str], n_groups: int = 2) -> list[list[str]]:
    """Split data rows into ``n_groups`` contiguous, roughly equal groups."""
    size = -(-len(rows) // n_groups)  # ceil division
    return [rows[i:i + size] for i in range(0, len(rows), size)]


def fdr_axis_top(df: pd.DataFrame, tfdr: float) -> float:
    """Upper FDR y-limit: the larger of the data max and the tFDR line, +headroom."""
    fdr_max = df.loc[df["metric"] == "FDR", "value"].max()
    return max(fdr_max, tfdr if tfdr > 0 else 0.0) * Y_HEADROOM


def choose_ticks(values: list[float], max_ticks: int = MAX_XTICKS) -> list[float]:
    """Return a readable subset of x-tick positions for a (log-scaled) SNR grid.

    Small grids (e.g. demo 02's 5 SNR levels) keep every value. Dense grids (e.g.
    demo 03's 21 levels) are thinned to ~``max_ticks`` values spaced evenly in
    log space -- matching the log x-axis -- so labels do not pile up at the high
    end where linearly-spaced values compress together. Endpoints are retained.
    """
    import math

    values = sorted(values)
    if len(values) <= max_ticks:
        return values

    lo, hi = math.log10(values[0]), math.log10(values[-1])
    ticks: list[float] = []
    for i in range(max_ticks):
        target = lo + (hi - lo) * i / (max_ticks - 1)
        nearest = min(values, key=lambda v: abs(math.log10(v) - target))
        if nearest not in ticks:
            ticks.append(nearest)
    return ticks


def marker_stride(values: list[float]) -> int:
    """Show at most ~10 markers per curve so dense grids stay legible."""
    return max(1, len(values) // 10)


def plot_fdr_tpr_vs_snr(
    df: pd.DataFrame,
    title: str,
    tfdr: float = DEFAULT_TFDR,
    legend_title: str = DEFAULT_LEGEND_TITLE,
) -> plt.Figure:
    """Pattern 1 — side-by-side FDR and TPR versus SNR panels.

    When ``tfdr > 0`` a black dotted reference line marks the target FDR level on
    the FDR panel, and the FDR y-axis is scaled to include it.
    """
    rows = list(dict.fromkeys(df["solver"]))
    snr_values = sorted(df["snr"].unique())
    colors = row_colors(rows)
    xticks = choose_ticks(snr_values)
    mev = marker_stride(snr_values)

    fig, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(15.0, 6.0),
        sharex=True,
        constrained_layout=False,
    )

    legend_handles = []
    legend_labels = []

    for ax, metric in zip(axes, METRICS):
        metric_df = df.loc[df["metric"] == metric]

        if metric == "FDR" and tfdr > 0:
            tfdr_line = ax.axhline(
                tfdr,
                color=TFDR_COLOR,
                linestyle=":",
                linewidth=2.0,
                zorder=1,
                label=f"tFDR = {tfdr:g}",
            )
            legend_handles.append(tfdr_line)
            legend_labels.append(f"tFDR = {tfdr:g}")

        for index, row in enumerate(rows):
            row_df = metric_df.loc[
                metric_df["solver"] == row
            ].sort_values("snr")

            if row_df.empty:
                continue

            line, = ax.plot(
                row_df["snr"],
                row_df["value"],
                color=colors[row],
                marker=MARKERS[index % len(MARKERS)],
                markersize=6.5,
                markevery=mev,
                markeredgecolor="white",
                markeredgewidth=0.7,
                linewidth=2.0,
                label=row,
            )

            if metric == "FDR":
                legend_handles.append(line)
                legend_labels.append(row)

        ax.set_xscale("log")
        ax.set_xticks(xticks)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_xlim(min(snr_values) * 0.82, max(snr_values) * 1.18)
        ax.set_xlabel("Signal-to-noise ratio (SNR)", fontsize=12)
        ax.set_ylabel(METRIC_LABELS[metric], fontsize=12)
        ax.set_title(metric, fontsize=14, fontweight="bold")
        ax.grid(True, which="major", linestyle="--", linewidth=0.7, alpha=0.45)
        ax.tick_params(axis="both", labelsize=11)

        if metric == "FDR":
            # Include the tFDR line in the range, plus a headroom budget.
            ax.set_ylim(0.0, fdr_axis_top(df, tfdr))
        else:
            ax.set_ylim(-0.02, 1.02)

    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.98)

    fig.legend(
        legend_handles,
        legend_labels,
        title=legend_title,
        loc="center left",
        bbox_to_anchor=(0.995, 0.5),
        frameon=True,
        fontsize=10,
        title_fontsize=11,
        ncol=1,
    )

    fig.subplots_adjust(
        left=0.08, right=0.78, bottom=0.12, top=0.86, wspace=0.28
    )
    return fig


def plot_fdr_tpr_vs_snr_grouped(
    df: pd.DataFrame,
    title: str,
    tfdr: float = DEFAULT_TFDR,
    legend_title: str = DEFAULT_LEGEND_TITLE,
    n_groups: int = 2,
) -> plt.Figure:
    """Pattern 2 — 2xN layout to de-clutter: rows = FDR/TPR, columns = row groups.

    Each panel carries only a fraction of the curves, so overlapping regions are
    easier to read than in the combined overview. Colours and markers are
    identical to the overview figure, and all FDR panels share the same y-scale
    for comparison.
    """
    rows = list(dict.fromkeys(df["solver"]))
    snr_values = sorted(df["snr"].unique())
    colors = row_colors(rows)
    markers = {s: MARKERS[i % len(MARKERS)] for i, s in enumerate(rows)}
    groups = split_rows(rows, n_groups=n_groups)
    xticks = choose_ticks(snr_values)
    mev = marker_stride(snr_values)

    # Row index ranges for the column titles (1-based, inclusive).
    bounds = []
    start = 1
    for group in groups:
        bounds.append((start, start + len(group) - 1))
        start += len(group)

    fig, axes = plt.subplots(
        nrows=2,
        ncols=len(groups),
        figsize=(7.5 * len(groups), 10.0),
        sharex=True,
        layout="constrained",
        squeeze=False,
    )

    for row_idx, metric in enumerate(METRICS):
        metric_df = df.loc[df["metric"] == metric]
        for col, group in enumerate(groups):
            ax = axes[row_idx, col]

            if metric == "FDR" and tfdr > 0:
                ax.axhline(
                    tfdr, color=TFDR_COLOR, linestyle=":", linewidth=2.0,
                    zorder=1, label=f"tFDR = {tfdr:g}",
                )

            for row in group:
                row_df = metric_df.loc[
                    metric_df["solver"] == row
                ].sort_values("snr")
                if row_df.empty:
                    continue
                ax.plot(
                    row_df["snr"], row_df["value"],
                    color=colors[row], marker=markers[row],
                    markersize=6.5, markevery=mev, markeredgecolor="white",
                    markeredgewidth=0.7, linewidth=2.0, label=row,
                )

            ax.set_xscale("log")
            ax.set_xticks(xticks)
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.set_xlim(min(snr_values) * 0.82, max(snr_values) * 1.18)
            ax.grid(
                True, which="major", linestyle="--", linewidth=0.7, alpha=0.45
            )
            ax.tick_params(axis="both", labelsize=11)
            ax.legend(
                fontsize=9,
                loc="upper left" if metric == "FDR" else "lower right",
                framealpha=0.9,
            )

            if metric == "FDR":
                ax.set_ylim(0.0, fdr_axis_top(df, tfdr))
            else:
                ax.set_ylim(-0.02, 1.02)

            if col == 0:
                ax.set_ylabel(METRIC_LABELS[metric], fontsize=12)
            if row_idx == 0:
                ax.set_title(
                    f"{legend_title} {bounds[col][0]}–{bounds[col][1]}",
                    fontsize=13, fontweight="bold",
                )
            if row_idx == len(METRICS) - 1:
                ax.set_xlabel("Signal-to-noise ratio (SNR)", fontsize=12)

    fig.suptitle(title, fontsize=15, fontweight="bold")
    return fig


def write_fdr_tpr_vs_snr_html(
    df: pd.DataFrame,
    title: str,
    out_html: Path,
    tfdr: float = DEFAULT_TFDR,
    legend_title: str = DEFAULT_LEGEND_TITLE,
    inline_js: bool = True,
) -> bool:
    """Pattern 3 — self-contained interactive HTML (FDR/TPR subplots, hover, toggles).

    Clicking a row in the legend hides/shows its curve in *both* panels, which
    is the interactive way to de-clutter. Returns False (with a note) if plotly is
    not installed.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        print(
            "note: plotly not installed - skipping interactive HTML "
            "(pip install plotly, or pass --no-plotly to silence)."
        )
        return False

    rows = list(dict.fromkeys(df["solver"]))
    snr_values = sorted(df["snr"].unique())
    colors = row_colors(rows)

    fig = make_subplots(
        rows=1, cols=2, subplot_titles=("FDR", "TPR"), horizontal_spacing=0.08
    )

    for metric, col in (("FDR", 1), ("TPR", 2)):
        metric_df = df.loc[df["metric"] == metric]
        for row in rows:
            row_df = metric_df.loc[
                metric_df["solver"] == row
            ].sort_values("snr")
            if row_df.empty:
                continue
            fig.add_trace(
                go.Scatter(
                    x=row_df["snr"],
                    y=row_df["value"],
                    mode="lines+markers",
                    name=row,
                    legendgroup=row,               # toggle both panels together
                    showlegend=(metric == "FDR"),  # list each row once
                    line=dict(color=rgba_to_css(colors[row]), width=2),
                    marker=dict(size=7),
                    hovertemplate=(
                        f"<b>{row}</b><br>SNR=%{{x}}<br>"
                        f"{metric}=%{{y:.4f}}<extra></extra>"
                    ),
                ),
                row=1,
                col=col,
            )

    if tfdr > 0:
        fig.add_hline(
            y=tfdr,
            line=dict(color="black", dash="dot", width=2),
            annotation_text=f"tFDR = {tfdr:g}",
            annotation_position="top left",
            row=1,
            col=1,
        )

    for col in (1, 2):
        fig.update_xaxes(
            type="log",
            tickvals=choose_ticks(snr_values),
            title_text="Signal-to-noise ratio (SNR)",
            row=1,
            col=col,
        )
    fig.update_yaxes(title_text="False discovery rate", row=1, col=1)
    fig.update_yaxes(
        title_text="True positive rate", range=[-0.02, 1.02], row=1, col=2
    )
    fig.update_layout(
        title=title.replace("$", ""),  # plain text: avoid raw LaTeX in HTML
        template="plotly_white",
        height=600,
        legend=dict(title=f"{legend_title} (click to toggle)"),
        hovermode="closest",
    )

    out_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(out_html, include_plotlyjs=(True if inline_js else "cdn"))
    print(f"Wrote: {out_html}")
    return True


def save_figure(
    fig: plt.Figure, outdir: Path, stem: str, formats: list[str]
) -> None:
    """Write the figure to outdir/<stem>.<ext> for each requested format."""
    outdir.mkdir(parents=True, exist_ok=True)
    for extension in formats:
        output_path = outdir / f"{stem}.{extension}"
        fig.savefig(
            output_path,
            dpi=300 if extension == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )
        print(f"Wrote: {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "csv_path",
        type=Path,
        help="Path to the tidy results CSV (solver,metric,snr,value).",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Figure super-title (default: derived from the CSV file name).",
    )
    parser.add_argument(
        "--legend-title",
        default=DEFAULT_LEGEND_TITLE,
        help=(
            "Legend title naming what the CSV's 'solver' column holds "
            f"(default: {DEFAULT_LEGEND_TITLE!r}; e.g. 'L-strategy' for demo "
            "04, 'Dummy distribution' for demo 05)."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help=(
            "Output directory (default: the simulation_results/plots/ sibling "
            "when the CSV lives in simulation_results/data/, else the CSV's "
            "own directory)."
        ),
    )
    parser.add_argument(
        "--stem",
        default=None,
        help=(
            "Output filename stem (default: <csv-stem>" + STEM_SUFFIX + ")."
        ),
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf"],
        help="Output formats to write (default: png pdf).",
    )
    parser.add_argument(
        "--tfdr",
        type=float,
        default=DEFAULT_TFDR,
        help=(
            "Target FDR level drawn as a dotted reference line on the FDR "
            f"panel (default: {DEFAULT_TFDR:g}). Pass 0 to omit the line."
        ),
    )
    parser.add_argument(
        "--groups",
        type=int,
        default=2,
        help="Number of row groups in the grouped figure (default: 2).",
    )
    parser.add_argument(
        "--no-grouped",
        action="store_true",
        help="Skip the grouped static figure (<stem>_grouped).",
    )
    parser.add_argument(
        "--no-plotly",
        action="store_true",
        help="Skip the interactive Plotly HTML (<stem>.html).",
    )
    parser.add_argument(
        "--plotly-cdn",
        action="store_true",
        help=(
            "Load plotly.js from a CDN instead of inlining it: a tiny HTML "
            "file that needs internet to render (default: self-contained)."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    csv_path = args.csv_path.resolve()

    df = read_mc_results(csv_path)
    title = (
        args.title
        if args.title is not None
        else f"T-Rex Monte Carlo simulation ({csv_path.stem})"
    )
    outdir = args.outdir if args.outdir is not None else default_plots_dir(csv_path)
    stem = args.stem if args.stem is not None else csv_path.stem + STEM_SUFFIX

    outdir = outdir.resolve()

    # 1. Combined overview (FDR | TPR).
    fig = plot_fdr_tpr_vs_snr(df, title, tfdr=args.tfdr,
                              legend_title=args.legend_title)
    save_figure(fig, outdir, stem, args.formats)
    plt.close(fig)

    # 2. De-cluttered grouped view.
    if not args.no_grouped:
        grouped_fig = plot_fdr_tpr_vs_snr_grouped(
            df, title, tfdr=args.tfdr,
            legend_title=args.legend_title, n_groups=args.groups,
        )
        save_figure(grouped_fig, outdir, stem + "_grouped", args.formats)
        plt.close(grouped_fig)

    # 3. Interactive HTML.
    if not args.no_plotly:
        write_fdr_tpr_vs_snr_html(
            df, title, outdir / f"{stem}.html",
            tfdr=args.tfdr, legend_title=args.legend_title,
            inline_js=not args.plotly_cdn,
        )


if __name__ == "__main__":
    main()
