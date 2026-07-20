#!/usr/bin/env python3
"""
trex_plt_utils.py — shared plotting utilities for the classical T-Rex demo suite.

Suite-level companion of trex_sim_utils.hpp: one module hosting the plot
patterns for every demo in this folder whose results use the tidy CSV schema

    solver,metric,snr,value

(the "solver" column holds whatever the demo sweeps: base solvers in demos
02/03, L-loop strategies in demo 04, dummy distributions in demo 05, storage
scenarios in demo 07, ...).
Row colours are inferred from the data, so no per-demo configuration is
required; each demo's generate_plots.sh supplies its CSV path and titles.

Plot patterns hosted here:
    plot_fdr_tpr_vs_snr()          — side-by-side FDR | TPR panels (static).
    plot_fdr_tpr_vs_snr_grouped()  — de-cluttered 2x2 grid, rows split in half.
    write_fdr_tpr_vs_snr_html()    — interactive Plotly HTML (hover, toggles).
    plot_comparison_grid()         — cross-solver 2xN grid: FDR/TPR rows x one
                                     column per base solver, every sweep row as a
                                     line (demo 04 L-strategies, demo 05 dummy
                                     distributions -- one column per TLARS/TOMP/
                                     TAFS on a shared scale).
    plot_scalability_dashboard()   — demo 08 runtime / memory / FDR / TPR vs
                                     problem size (its own wider CSV schema,
                                     read via read_scalability()).

Three CLI modes (all typically driven from a demo's generate_plots.sh):
    # default -- per-CSV overview + grouped + interactive html
    python trex_plt_utils.py <csv> --title "..." --legend-title "L-strategy"
    # cross-solver comparison grid from several per-solver CSVs
    python trex_plt_utils.py grid --labels TLARS TOMP TAFS \
        --csvs a_tlars.csv a_tomp.csv a_tafs.csv \
        --legend-title "L-strategy" --stem my_grid --outdir <plots>
    # demo-08 scalability dashboard from the wide CSV
    python trex_plt_utils.py scalability <wide_csv> --stem my_dash --outdir <plots>

Output location: figures default to the simulation_results/plots/ sibling when
the CSV lives in simulation_results/data/ (the suite convention); otherwise
they are written next to the CSV. Override with --outdir. (The grid and
scalability modes require an explicit --outdir, since they are not tied to a
single CSV's location.)

All three modes accept ``--xscale {linear,sqrt,log}`` (default ``linear``,
matching the rest of the suites): ``sqrt`` spreads a sweep grid that is dense
near zero while keeping 0 on-axis, ``log`` restores the logarithmic axis with
the log-spaced tick values from choose_ticks().
"""

import argparse
import re
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any, cast

import matplotlib

matplotlib.use("Agg")  # headless / batch-safe: no display required

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, ScalarFormatter

# ---------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------
STEM_SUFFIX = "_fdr_tpr_vs_snr"
DEFAULT_TFDR = 0.10          # target FDR control level (tFDR)
# Names what the CSV's 'solver' column holds. "Method" is the suite-wide default
# (matching trex_da_plt_utils.py), since for most demos the column lists the
# selector methods themselves. Demos whose column holds something else pass their
# own: 'L-strategy' (demo 04), 'Dummy distribution' (demo 05), 'Scenario'
# (demo 07's storage overlay).
DEFAULT_LEGEND_TITLE = "Method"
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


def natural_sort_key(name: str) -> list:
    """Sort key that splits digit runs from text so numbered labels order
    numerically: ``SKIPL_5p`` < ``SKIPL_10p`` < ``SKIPL_20p`` (a plain string
    sort would give 10p, 20p, 5p). Non-numeric parts compare case-insensitively.
    """
    return [
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", str(name))
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

    # Order rows by a *natural* sort of the 'solver' column so numbered rows
    # read numerically -- e.g. demo 04's SKIPL_5p < SKIPL_10p < SKIPL_20p, where
    # a plain string sort would give 10p, 20p, 5p. This ordering flows straight
    # into the legend and the colour assignment. For row sets without embedded
    # numbers it matches the previous plain alphabetical order.
    rank = {
        name: index
        for index, name in enumerate(
            sorted(df["solver"].unique(), key=natural_sort_key)
        )
    }
    df = df.assign(_row_rank=df["solver"].map(rank))
    return (
        df.sort_values(["metric", "_row_rank", "snr"])
        .drop(columns="_row_rank")
    )


def default_plots_dir(csv_path: Path) -> Path:
    """Suite convention: CSVs live in simulation_results/data/, figures in the
    sibling simulation_results/plots/. Falls back to the CSV's own directory
    for CSVs outside that layout."""
    if csv_path.parent.name == "data":
        return csv_path.parent.parent / "plots"
    return csv_path.parent


def qualitative_colors(name: str) -> list:
    """The hand-picked colour list behind a qualitative colormap.

    ``get_cmap`` is typed as returning the base ``Colormap``, but ``.colors``
    only exists on ``ListedColormap`` -- which is what the qualitative maps
    (``tab10``, ``Dark2``) actually are. Narrowing here keeps the read type-safe
    and fails loudly on a colormap that carries no discrete colour list at all
    (a ``LinearSegmentedColormap`` such as jet).
    """
    cmap = plt.get_cmap(name)
    if not isinstance(cmap, ListedColormap):
        raise TypeError(
            f"{name!r} is a {type(cmap).__name__}, which has no discrete "
            "colour list; pass a ListedColormap such as tab10 or Dark2."
        )
    # `.colors` is declared as the broad ArrayLike (which admits scalars); on a
    # ListedColormap it is concretely the sequence of RGB tuples.
    return list(cast(Sequence[Any], cmap.colors))


def row_colors(rows: list[str]) -> dict[str, tuple]:
    """Assign a stable, clearly distinct colour to each data row.

    Concatenates two *qualitative* colormaps (``tab10`` + ``Dark2`` = 18 hand-
    picked, high-contrast colours) rather than sampling a smooth gradient, so
    adjacent curves are easy to tell apart. Colours are keyed by the row's
    position in the (file-order) list, so the mapping is deterministic and works
    for any demo's row set without a hardcoded table.
    """
    palette = qualitative_colors("tab10") + qualitative_colors("Dark2")
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


DEFAULT_XSCALE = "linear"
# The demo-08 dashboard's x-axis is a problem size p spanning decades, not a
# sweep grid, so it keeps a logarithmic default (see _scal_p_axis).
SCAL_XSCALE = "log"
XSCALES = ("linear", "sqrt", "log")


def apply_sweep_xaxis(ax, values: list[float], xscale: str = DEFAULT_XSCALE,
                      ticks: list[float] | None = None,
                      force_ticks: bool = False) -> None:
    """Configure a sweep axis' scale, ticks and limits.

    ``linear`` (the default, matching the DA-TRex / Screen-TRex / T-Rex+GVS
    suites) and ``sqrt`` take their major ticks from the locators rather than
    from the swept values, so labels cannot collide however tightly the grid is
    spaced, and both carry minor ticks with a faint minor x-grid. ``sqrt``
    stretches the region near zero while keeping 0 on-axis. ``log`` is the
    logarithmic axis this module used to hard-code, tick-marked at ``ticks``
    (the log-spaced subset from choose_ticks()).

    ``ticks`` is honoured on the log axis, and on every scale when
    ``force_ticks`` is set -- used by the scalability p-axis, whose tick *labels*
    carry the paired sample size and so must line up with fixed positions.

    Mirrors ``_apply_sweep_xaxis`` of ``../trex_da/trex_da_plt_utils.py``.
    """
    lo, hi = min(values), max(values)
    if xscale == "log":
        ax.set_xscale("log")
        if ticks is not None:
            ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.set_xlim(lo * 0.82, hi * 1.18)
        return

    if xscale == "sqrt":
        ax.set_xscale(
            "function",
            functions=(
                lambda x: np.sqrt(np.maximum(x, 0.0)),
                lambda x: np.square(x),
            ),
        )
        span_sqrt = np.sqrt(hi) - np.sqrt(max(lo, 0.0))
        pad = (span_sqrt * 0.04) ** 2 if span_sqrt > 0 else 0.5
        ax.set_xlim(lo, hi + pad)
    else:
        span = hi - lo
        pad = span * 0.04 if span > 0 else 0.5
        ax.set_xlim(lo - pad, hi + pad)

    if force_ticks and ticks is not None:
        ax.set_xticks(ticks)
    else:
        ax.xaxis.set_major_locator(
            MaxNLocator(nbins="auto", steps=[1, 2, 2.5, 5, 10])
        )
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.grid(which="minor", axis="x", alpha=0.25, linewidth=0.5)


def plot_fdr_tpr_vs_snr(
    df: pd.DataFrame,
    title: str,
    tfdr: float = DEFAULT_TFDR,
    legend_title: str = DEFAULT_LEGEND_TITLE,
    xscale: str = DEFAULT_XSCALE,
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

        apply_sweep_xaxis(ax, snr_values, xscale, ticks=xticks)
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

    # This overview carries every row at once, so with many solvers the legend
    # would obscure the curves if placed inside a panel; it is kept just outside
    # the right (TPR) panel instead. Anchoring it a hair past the panel's right
    # edge (right=0.80 -> anchor 0.81) keeps it hugging the figure rather than
    # drifting off to the far right, and an opaque white box (framealpha=1.0)
    # avoids any see-through legend background.
    legend = fig.legend(
        legend_handles,
        legend_labels,
        title=legend_title,
        loc="center left",
        bbox_to_anchor=(0.81, 0.5),
        frameon=True,
        framealpha=1.0,
        facecolor="white",
        edgecolor="0.7",
        fontsize=10,
        title_fontsize=11,
        ncol=1,
    )
    legend.get_frame().set_linewidth(0.8)

    fig.subplots_adjust(
        left=0.08, right=0.80, bottom=0.12, top=0.86, wspace=0.28
    )
    return fig


def plot_fdr_tpr_vs_snr_grouped(
    df: pd.DataFrame,
    title: str,
    tfdr: float = DEFAULT_TFDR,
    legend_title: str = DEFAULT_LEGEND_TITLE,
    n_groups: int = 2,
    xscale: str = DEFAULT_XSCALE,
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

            apply_sweep_xaxis(ax, snr_values, xscale, ticks=xticks)
            ax.grid(
                True, which="major", linestyle="--", linewidth=0.7, alpha=0.45
            )
            ax.tick_params(axis="both", labelsize=11)
            leg = ax.legend(
                fontsize=9,
                loc="upper left" if metric == "FDR" else "lower right",
                frameon=True,
                framealpha=1.0,
                facecolor="white",
                edgecolor="0.7",
            )
            leg.get_frame().set_linewidth(0.8)

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
    xscale: str = DEFAULT_XSCALE,
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

    # The html and the png/pdf of the same figure share an x-axis. Plotly has no
    # sqrt axis, so it falls back to a linear one with minor ticks; the ticks
    # come from plotly's autorange there, so near-neighbours cannot overlap.
    if xscale == "log":
        xaxis_kwargs: dict = dict(type="log", tickvals=choose_ticks(snr_values))
    else:
        xaxis_kwargs = dict(type="linear", minor=dict(ticklen=4, showgrid=True))
    for col in (1, 2):
        fig.update_xaxes(
            title_text="Signal-to-noise ratio (SNR)",
            row=1,
            col=col,
            **xaxis_kwargs,
        )
    fig.update_yaxes(title_text="False discovery rate", row=1, col=1)
    fig.update_yaxes(
        title_text="True positive rate", range=[-0.02, 1.02], row=1, col=2
    )

    # Draw a full frame (box) around every panel. mirror=True repeats each axis
    # line on the opposite side of its own panel, so the bottom/top (x) and
    # left/right (y) lines close the box -- matching the framed look of the
    # static matplotlib panels (plotly_white otherwise leaves top/right open).
    # Applied globally so it covers both subplots without clobbering the
    # per-axis titles/ranges set above.
    fig.update_xaxes(
        showline=True, linewidth=1, linecolor="rgb(82,82,82)", mirror=True
    )
    fig.update_yaxes(
        showline=True, linewidth=1, linecolor="rgb(82,82,82)", mirror=True
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


def ordered_grid_rows(frames: dict) -> list[str]:
    """Union of the 'solver'-column values across per-column frames, in
    first-seen order. read_mc_results already natural-sorts each frame's rows
    (so e.g. SKIPL_5p < SKIPL_10p < SKIPL_20p); keeping one ordered list makes
    the colour / marker assignment identical in every column, which is the whole
    point of the cross-solver view.
    """
    seen: list[str] = []
    for df in frames.values():
        for name in df["solver"]:  # the swept row: L-strategy / distribution
            if name not in seen:
                seen.append(name)
    return seen


def plot_comparison_grid(
    frames: dict,
    title: str,
    legend_title: str = DEFAULT_LEGEND_TITLE,
    tfdr: float = DEFAULT_TFDR,
    xscale: str = DEFAULT_XSCALE,
) -> plt.Figure:
    """Pattern 4 — cross-solver comparison grid.

    ``frames`` is an *ordered* mapping ``{column_label: tidy_df}`` (one per base
    solver, e.g. TLARS / TOMP / TAFS). The figure is a 2xN grid: rows = metric
    (FDR on top, TPR below), columns = the frames' keys. Every panel carries all
    of the swept rows (the CSV's "solver" column -- L-strategies in demo 04,
    dummy distributions in demo 05) as lines vs SNR, with colours and markers
    keyed to the row so a given row is the same colour in every column. All FDR
    panels share one y-scale and all TPR panels share [0, 1], so the columns are
    directly comparable. The (potentially long) row legend is drawn once, just
    right of the grid, with an opaque box.
    """
    columns = list(frames.keys())
    rows = ordered_grid_rows(frames)
    colors = row_colors(rows)
    markers = {r: MARKERS[i % len(MARKERS)] for i, r in enumerate(rows)}

    # A shared SNR grid / tick set / y-top across panels keeps them comparable.
    all_snr = sorted(
        {snr for df in frames.values() for snr in df["snr"].unique()}
    )
    xticks = choose_ticks(all_snr)
    mev = marker_stride(all_snr)
    fdr_top = max(fdr_axis_top(df, tfdr) for df in frames.values())

    fig, axes = plt.subplots(
        nrows=2,
        ncols=len(columns),
        figsize=(5.4 * len(columns) + 2.4, 9.0),
        sharex=True,
        squeeze=False,
    )

    handles: list = []
    labels: list = []

    for row_idx, metric in enumerate(METRICS):
        for col, column in enumerate(columns):
            ax = axes[row_idx, col]
            metric_df = frames[column].loc[frames[column]["metric"] == metric]

            if metric == "FDR" and tfdr > 0:
                tfdr_line = ax.axhline(
                    tfdr, color=TFDR_COLOR, linestyle=":", linewidth=2.0,
                    zorder=1, label=f"tFDR = {tfdr:g}",
                )
                if not handles:  # capture the reference line once, first
                    handles.append(tfdr_line)
                    labels.append(f"tFDR = {tfdr:g}")

            for row in rows:
                row_df = metric_df.loc[
                    metric_df["solver"] == row
                ].sort_values("snr")
                if row_df.empty:
                    continue
                line, = ax.plot(
                    row_df["snr"], row_df["value"],
                    color=colors[row], marker=markers[row],
                    markersize=6.0, markevery=mev, markeredgecolor="white",
                    markeredgewidth=0.7, linewidth=1.9, label=row,
                )
                # Build the shared legend once, from the top-left panel.
                if row_idx == 0 and col == 0:
                    handles.append(line)
                    labels.append(row)

            apply_sweep_xaxis(ax, all_snr, xscale, ticks=xticks)
            ax.grid(
                True, which="major", linestyle="--", linewidth=0.7, alpha=0.45
            )
            ax.tick_params(axis="both", labelsize=11)

            if metric == "FDR":
                ax.set_ylim(0.0, fdr_top)
            else:
                ax.set_ylim(-0.02, 1.02)

            if col == 0:
                ax.set_ylabel(METRIC_LABELS[metric], fontsize=12)
            if row_idx == 0:
                ax.set_title(column, fontsize=14, fontweight="bold")
            if row_idx == len(METRICS) - 1:
                ax.set_xlabel("Signal-to-noise ratio (SNR)", fontsize=12)

    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.99)

    legend = fig.legend(
        handles, labels,
        title=legend_title,
        loc="center left",
        bbox_to_anchor=(0.845, 0.5),
        frameon=True, framealpha=1.0, facecolor="white", edgecolor="0.7",
        fontsize=10, title_fontsize=11, ncol=1,
    )
    legend.get_frame().set_linewidth(0.8)

    fig.subplots_adjust(
        left=0.07, right=0.84, bottom=0.08, top=0.90, wspace=0.18, hspace=0.12
    )
    return fig


# ---------------------------------------------------------------------
# Demo 08 — scalability dashboard (wider CSV schema than the rest)
# ---------------------------------------------------------------------
# Demo 08 is a 2-D sweep (scenario x solver over the SNR grid) on the fully
# memory-mapped pipeline; its CSV uses a wider schema than the tidy suite one,
#     scenario,solver,n,p,num_mc,snr,metric,value
# so it is read via read_scalability() (not read_mc_results) and rendered by
# plot_scalability_dashboard().

# Base solvers, in a fixed display order with stable colours / markers. These
# keys are the raw names the demo binary writes into the CSV's 'solver' column.
SCAL_SOLVERS = ("TLARS", "TOMP", "TAFS")
SCAL_SOLVER_COLORS = {"TLARS": "#1f77b4", "TOMP": "#d62728", "TAFS": "#2ca02c"}
SCAL_SOLVER_MARKERS = {"TLARS": "o", "TOMP": "s", "TAFS": "^"}

# TAFS is the only solver here run with a non-zero AFS correlation parameter;
# the value mirrors the solver table in demo_trex_08_mc_sim_scalability.cpp
# ({TAFS, "TAFS", lambda2=0.0, rho_afs=0.3}) -- keep the two in sync. It is a
# defining part of the method rather than a tuning detail, so it is named in the
# legend. Applied to the *label* only: the lookups above stay keyed on the raw
# CSV name.
SCAL_SOLVER_RHO_AFS = {"TAFS": 0.3}


def scal_solver_label(solver: str) -> str:
    """Legend label for a scalability solver, annotated with its defining
    hyper-parameter where it has one (``TAFS`` -> ``TAFS (rho = 0.3)``).

    Plain ASCII rather than mathtext, matching the DA suite's labels.
    """
    rho = SCAL_SOLVER_RHO_AFS.get(solver)
    return solver if rho is None else f"{solver} (rho = {rho:g})"

# Column labels as they appear in the CSV's 'metric' field.
M_RUNTIME = "Total s"
M_RSS = "RSS MiB"        # resident set size after select() (MiB) - tracks data
M_PEAKRSS = "PeakRSS"    # process peak RSS (MiB) - inflated by cold-heap footprint
M_XY = "X+y MiB"         # on-disk X/y mmap backing size (MiB)
M_FDR = "FDR"
M_TPR = "TPR"
SCAL_REF_COLOR = "0.35"  # X/y backing-size reference line


def read_scalability(path: Path) -> pd.DataFrame:
    """Load the wide demo-08 CSV, normalising column and metric whitespace."""
    if not path.exists():
        raise FileNotFoundError(
            f"Input file not found:\n  {path}\n"
            "Pass the results CSV path as the first argument."
        )
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip().str.lower()

    required = {"scenario", "solver", "n", "p", "snr", "metric", "value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"CSV is missing required columns: {sorted(missing)}\n"
            f"Found columns: {df.columns.tolist()}"
        )

    df["solver"] = df["solver"].astype(str).str.strip()
    df["metric"] = df["metric"].astype(str).str.strip()
    for col in ("n", "p", "snr", "value"):
        df[col] = pd.to_numeric(df[col], errors="raise")
    return df


def _scal_metric_by_p(df: pd.DataFrame, solver: str, metric: str):
    """Mean / min / max of a metric over the SNR grid, per problem size p.

    Returns (p_values, mean, lo, hi) sorted by p. Aggregating over SNR collapses
    the third sweep axis into a band, so each solver is a single curve in p.
    """
    sub = df[(df["solver"] == solver) & (df["metric"] == metric)]
    if sub.empty:
        return [], [], [], []
    stats = sub.groupby("p")["value"].agg(["mean", "min", "max"]).sort_index()
    return (
        stats.index.tolist(), stats["mean"].tolist(),
        stats["min"].tolist(), stats["max"].tolist(),
    )


def _scal_metric_by_p_all(df: pd.DataFrame, metric: str):
    """Mean of a metric per problem size p, aggregated over *all* solvers and
    SNR levels. Used for the memory quantities, which are data-dominated and so
    essentially solver-independent -- one curve per quantity keeps the panel
    readable instead of three overlapping ones.
    """
    sub = df[df["metric"] == metric]
    if sub.empty:
        return [], []
    stats = sub.groupby("p")["value"].mean().sort_index()
    return stats.index.tolist(), stats.tolist()


def _scal_scatter_snr(ax, df: pd.DataFrame, solver: str, metric: str) -> None:
    """Overlay the individual per-SNR values as small translucent dots, so the
    spread behind a mean curve is visible without a dominating filled band."""
    sub = df[(df["solver"] == solver) & (df["metric"] == metric)]
    if not sub.empty:
        ax.scatter(
            sub["p"], sub["value"], s=18, color=SCAL_SOLVER_COLORS[solver],
            alpha=0.35, edgecolors="none", zorder=2,
        )


def _scal_p_axis(ax, p_values: list, n_for_p: dict,
                 xscale: str = SCAL_XSCALE) -> None:
    """x-axis in p, tick labels annotated with the paired sample size n.

    The ticks stay pinned to the measured p values on every scale (their labels
    carry the paired n, so a locator-chosen position would be mislabelled).
    Unlike the SNR sweeps this axis defaults to ``log``: p is a problem *size*
    spanning decades (1e3 ... 1e5), so on a linear axis the low decades collapse
    into each other and their two-line tick labels overlap. ``--xscale`` still
    overrides it.
    """
    apply_sweep_xaxis(ax, p_values, xscale, ticks=p_values, force_ticks=True)
    if xscale == "log":
        # Roomier than the shared log limits: the p grid is sparse.
        ax.set_xlim(min(p_values) / 1.6, max(p_values) * 1.6)
    ax.set_xticklabels(
        [f"{int(p):,}\n(n={n_for_p[int(p)]:,})" for p in p_values]
    )
    ax.set_xlabel("Number of features $p$  (paired sample size $n$)", fontsize=11)


def plot_scalability_dashboard(
    df: pd.DataFrame, title: str, tfdr: float = DEFAULT_TFDR,
    xscale: str = SCAL_XSCALE,
) -> plt.Figure:
    """Pattern 5 — demo-08 runtime / memory / quality scaling dashboard.

    A 2x2 grid, x-axis = problem size p (log), one line per base solver,
    aggregated over the SNR grid: runtime (min-max band), memory (backing size
    vs resident RSS vs process peak), and FDR / TPR at scale (mean line with
    translucent per-SNR dots).
    """
    solvers = [s for s in SCAL_SOLVERS if s in set(df["solver"])]
    p_values = sorted(df["p"].unique())
    n_for_p = {
        int(p): int(df.loc[df["p"] == p, "n"].iloc[0]) for p in p_values
    }

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14.0, 10.0))
    ax_rt, ax_mem = axes[0]
    ax_fdr, ax_tpr = axes[1]

    # --- Runtime (log-log), mean over SNR with a min-max band ---------------
    for solver in solvers:
        p, mean, lo, hi = _scal_metric_by_p(df, solver, M_RUNTIME)
        if not p:
            continue
        ax_rt.fill_between(p, lo, hi, color=SCAL_SOLVER_COLORS[solver], alpha=0.15)
        ax_rt.plot(
            p, mean, color=SCAL_SOLVER_COLORS[solver],
            marker=SCAL_SOLVER_MARKERS[solver], markersize=8,
            markeredgecolor="white", markeredgewidth=0.8, linewidth=2.2,
            label=scal_solver_label(solver),
        )
    ax_rt.set_yscale("log")
    ax_rt.set_ylabel("Total wall-clock per trial [s]", fontsize=12)
    ax_rt.set_title("Runtime scaling", fontsize=14, fontweight="bold")

    # --- Memory (log-log): backing size vs resident RSS vs process peak ------
    # These three are data-dominated (nearly solver-independent), so each is a
    # single curve aggregated over solvers. The story lives in the *gaps*:
    # resident RSS after select() hugs the on-disk X/y size (memory-mapping
    # works), while the process peak RSS runs ~10x higher -- the cold per-
    # experiment coefficient-path heap the demo flags as the scalability limit.
    p_xy, mean_xy = _scal_metric_by_p_all(df, M_XY)
    if p_xy:
        ax_mem.plot(
            p_xy, mean_xy, color=SCAL_REF_COLOR, linestyle="--", linewidth=2.0,
            marker="D", markersize=7, markerfacecolor="white",
            label="$X$/$y$ mmap backing (on disk)",
        )
    p_rss, mean_rss = _scal_metric_by_p_all(df, M_RSS)
    if p_rss:
        ax_mem.plot(
            p_rss, mean_rss, color="#1f77b4", linewidth=2.4, marker="o",
            markersize=8, markeredgecolor="white", markeredgewidth=0.8,
            label="Resident RSS after select()",
        )
    p_peak, mean_peak = _scal_metric_by_p_all(df, M_PEAKRSS)
    if p_peak:
        ax_mem.plot(
            p_peak, mean_peak, color="#d62728", linewidth=2.4, marker="^",
            markersize=8, markeredgecolor="white", markeredgewidth=0.8,
            label="Process peak RSS (cold-heap inflated)",
        )
    ax_mem.set_yscale("log")
    ax_mem.set_ylabel("Memory [MiB]", fontsize=12)
    ax_mem.set_title(
        "Memory scaling: mmap keeps resident RSS near the data size",
        fontsize=13, fontweight="bold",
    )

    # --- FDR at scale -------------------------------------------------------
    if tfdr > 0:
        ax_fdr.axhline(
            tfdr, color="black", linestyle=":", linewidth=2.0,
            label=f"tFDR = {tfdr:g}",
        )
    fdr_top = tfdr * 1.25 if tfdr > 0 else 0.1
    for solver in solvers:
        p, mean, _lo, hi = _scal_metric_by_p(df, solver, M_FDR)
        if not p:
            continue
        fdr_top = max(fdr_top, max(hi) * 1.15)
        _scal_scatter_snr(ax_fdr, df, solver, M_FDR)
        ax_fdr.plot(
            p, mean, color=SCAL_SOLVER_COLORS[solver],
            marker=SCAL_SOLVER_MARKERS[solver], markersize=8,
            markeredgecolor="white", markeredgewidth=0.8, linewidth=2.2,
            label=scal_solver_label(solver),
        )
    ax_fdr.set_ylim(0.0, fdr_top)
    ax_fdr.set_ylabel("False discovery rate", fontsize=12)
    ax_fdr.set_title("FDR at scale (mean line, per-SNR dots)",
                     fontsize=14, fontweight="bold")

    # --- TPR at scale -------------------------------------------------------
    # TPR varies strongly with SNR, so per-SNR dots (not a filled band, which
    # would swamp the panel at the low-SNR baseline scenario) show the spread.
    for solver in solvers:
        p, mean, _lo, _hi = _scal_metric_by_p(df, solver, M_TPR)
        if not p:
            continue
        _scal_scatter_snr(ax_tpr, df, solver, M_TPR)
        ax_tpr.plot(
            p, mean, color=SCAL_SOLVER_COLORS[solver],
            marker=SCAL_SOLVER_MARKERS[solver], markersize=8,
            markeredgecolor="white", markeredgewidth=0.8, linewidth=2.2,
            label=scal_solver_label(solver),
        )
    ax_tpr.set_ylim(-0.02, 1.05)
    ax_tpr.set_ylabel("True positive rate", fontsize=12)
    ax_tpr.set_title("TPR at scale (mean line, per-SNR dots)",
                     fontsize=14, fontweight="bold")

    for ax in (ax_rt, ax_mem, ax_fdr, ax_tpr):
        _scal_p_axis(ax, p_values, n_for_p, xscale)
        ax.grid(True, which="both", linestyle="--", linewidth=0.7, alpha=0.4)
        ax.tick_params(axis="both", labelsize=10)
        leg = ax.legend(
            fontsize=9, loc="best",
            frameon=True, framealpha=1.0, facecolor="white", edgecolor="0.7",
        )
        leg.get_frame().set_linewidth(0.8)

    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.99)
    fig.subplots_adjust(
        left=0.075, right=0.975, bottom=0.085, top=0.90, wspace=0.20, hspace=0.32
    )
    return fig


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
            "04, 'Dummy distribution' for demo 05, 'Scenario' for demo 07)."
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
        "--xscale",
        choices=XSCALES,
        default=DEFAULT_XSCALE,
        help=(
            f"X-axis scale (default: {DEFAULT_XSCALE!r}). 'sqrt' spreads a "
            "sweep grid that is dense near zero; 'log' restores the "
            "logarithmic SNR axis."
        ),
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


def main_fdr_tpr() -> None:
    """Default mode: per-CSV overview + grouped + interactive HTML."""
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
                              legend_title=args.legend_title,
                              xscale=args.xscale)
    save_figure(fig, outdir, stem, args.formats)
    plt.close(fig)

    # 2. De-cluttered grouped view.
    if not args.no_grouped:
        grouped_fig = plot_fdr_tpr_vs_snr_grouped(
            df, title, tfdr=args.tfdr,
            legend_title=args.legend_title, n_groups=args.groups,
            xscale=args.xscale,
        )
        save_figure(grouped_fig, outdir, stem + "_grouped", args.formats)
        plt.close(grouped_fig)

    # 3. Interactive HTML.
    if not args.no_plotly:
        write_fdr_tpr_vs_snr_html(
            df, title, outdir / f"{stem}.html",
            tfdr=args.tfdr, legend_title=args.legend_title,
            inline_js=not args.plotly_cdn, xscale=args.xscale,
        )


def main_grid(argv: list) -> None:
    """`grid` mode: cross-solver comparison grid from several per-solver CSVs.

    ``--labels`` and ``--csvs`` are parallel lists (one column per pair, in the
    given order); a listed CSV that does not exist yet is skipped with a note, so
    the grid renders whatever solvers have been run so far.
    """
    parser = argparse.ArgumentParser(
        prog="trex_plt_utils.py grid",
        description="Cross-solver FDR/TPR comparison grid.",
    )
    parser.add_argument(
        "--labels", nargs="+", required=True,
        help="Column labels, one per CSV (e.g. TLARS TOMP TAFS).",
    )
    parser.add_argument(
        "--csvs", nargs="+", required=True, type=Path,
        help="Per-solver tidy CSVs, parallel to --labels.",
    )
    parser.add_argument("--title", default="T-Rex cross-solver comparison")
    parser.add_argument("--legend-title", default=DEFAULT_LEGEND_TITLE)
    parser.add_argument("--stem", required=True, help="Output filename stem.")
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Output directory (the grid is not tied to one CSV's location).",
    )
    parser.add_argument("--formats", nargs="+", default=["png", "pdf"])
    parser.add_argument("--tfdr", type=float, default=DEFAULT_TFDR)
    parser.add_argument("--xscale", choices=XSCALES, default=DEFAULT_XSCALE)
    # generate_plots.sh forwards its "$@" to every plotter; ignore flags only
    # the default mode understands (--no-plotly, --no-grouped, ...).
    args, _ignored = parser.parse_known_args(argv)

    if len(args.labels) != len(args.csvs):
        raise SystemExit(
            f"--labels ({len(args.labels)}) and --csvs ({len(args.csvs)}) "
            "must have the same length."
        )

    frames = {}
    for label, csv in zip(args.labels, args.csvs):
        csv = csv.resolve()
        if not csv.exists():
            print(f"note: skipping {label} - CSV not found:\n  {csv}")
            continue
        frames[label] = read_mc_results(csv)

    if not frames:
        raise SystemExit(
            "No input CSVs found - run the demo first, then re-run this plotter."
        )

    fig = plot_comparison_grid(
        frames, args.title, legend_title=args.legend_title, tfdr=args.tfdr,
        xscale=args.xscale,
    )
    save_figure(fig, args.outdir.resolve(), args.stem, args.formats)
    plt.close(fig)


def main_scalability(argv: list) -> None:
    """`scalability` mode: demo-08 runtime / memory / quality dashboard."""
    parser = argparse.ArgumentParser(
        prog="trex_plt_utils.py scalability",
        description="Demo-08 scalability dashboard (wide CSV schema).",
    )
    parser.add_argument(
        "csv_path", type=Path,
        help="Wide demo-08 CSV (scenario,solver,n,p,num_mc,snr,metric,value).",
    )
    parser.add_argument(
        "--title",
        default="T-Rex scalability on the memory-mapped pipeline",
    )
    parser.add_argument("--stem", required=True, help="Output filename stem.")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--formats", nargs="+", default=["png", "pdf"])
    parser.add_argument("--tfdr", type=float, default=DEFAULT_TFDR)
    parser.add_argument(
        "--xscale", choices=XSCALES, default=SCAL_XSCALE,
        help=(f"Scale of the problem-size (p) x-axis (default: "
              f"{SCAL_XSCALE!r} -- p spans decades)."),
    )
    args, _ignored = parser.parse_known_args(argv)

    df = read_scalability(args.csv_path.resolve())
    fig = plot_scalability_dashboard(df, args.title, tfdr=args.tfdr,
                                     xscale=args.xscale)
    save_figure(fig, args.outdir.resolve(), args.stem, args.formats)
    plt.close(fig)


def main() -> None:
    """Dispatch on the first token: `grid` / `scalability` sub-modes, else the
    default per-CSV FDR/TPR mode (so bare `trex_plt_utils.py <csv> ...` still
    works unchanged for demos 02/03 and the per-solver calls of 04/05)."""
    argv = sys.argv[1:]
    if argv and argv[0] == "grid":
        main_grid(argv[1:])
    elif argv and argv[0] == "scalability":
        main_scalability(argv[1:])
    else:
        main_fdr_tpr()


if __name__ == "__main__":
    main()
