# ==============================================================================
# demo_trex_gvs_08.py
# ==============================================================================
#
# T-Rex+GVS block-equicorrelated benchmark: HAC-discovered vs oracle groups.
#
# A 100-block design (G = 100 blocks of block_size = 5 tiling p = 500 columns),
# with 10 active blocks. Four method variants are compared over a (rho x snr)
# grid, for each of three within-block truth scenarios:
#
#   M1 = EN  + HAC-clustered groups   (groups discovered from the data)
#   M2 = EN  + oracle groups          (true block labels supplied)
#   M3 = IEN + HAC-clustered groups
#   M4 = IEN + oracle groups
#
#   Truth scenarios: INDIVIDUAL (1 active var/block), REPRESENTATIVE (2-3/block),
#                    WHOLE (all 5/block).
#
# The central comparisons are M1-vs-M2 and M3-vs-M4: how much accuracy does
# automatic group discovery (HAC) cost relative to knowing the true groups?
#
# Metrics: coordinate TPR/FDR, block_FDR, a scenario-specific block metric
# (blk_hit for INDIVIDUAL/REPRESENTATIVE, full_blk for WHOLE), T_stop, purity,
# and M_found.
#
# GROUPING DIAGNOSTICS: the GVSSelectionResult returned by select() exposes
# `max_clusters` (number of clusters the selector formed -> M_found) and
# `groups_vec` (per-variable cluster labels). Both `purity` (the block-purity
# rate, per the cpp benchmark) and `M_found` are computed directly from that
# result, for HAC and oracle methods alike. This matches the cpp benchmark and
# surfaces more than the R demo: the R6 API exposes only `$max_clusters`, so the
# R demo reports purity as 1.0/"-" rather than the real per-block rate.
#
# Console-only (mirrors cpp/R demo 08). Part 1 (single-run) runs by default; the
# full benchmark (Part 2) is opt-in via RUN_PART_2.
# Indices are 0-based (Python convention); the R counterpart is 1-based.
# ==============================================================================

import sys
import os
import multiprocessing

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

from trex_gvs_dgps import dgp_block_equicorr
from gvs_sim_common import N_CORES, build_gvs_selector
from trex_selector_neo.utils import compute_fdp, compute_tpp

N_WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else N_CORES

CFG = dict(
    n=200, p=500,
    G=100, block_size=5, n_active_blocks=10,
    tFDR=0.1, K=20,
    corr_max=0.5,           # much lower than the 0.98 used in demos 01-07
    b=3.0,
    base_seed=2026, num_MC=500,
    en_solver="TENET_AUG",
    snr_grid=[0.1, 0.2, 0.5, 1.0, 2.0, 5.0],
    rho_grid=[0.5, 0.9],
)

RUN_PART_1 = True    # single-run demo
RUN_PART_2 = False   # full benchmark (3 scenarios x rho x snr grid)

_METRIC_NAMES = ("coord_tpr", "coord_fdp", "block_hit", "block_fdp",
                 "full_block", "purity", "T_stop", "M_found")


# ==============================================================================
# Block-level metric helpers  (mirror the cpp/R benchmark; 0-based)
# ==============================================================================

def _block_cols(blk, block_size):
    """0-based column indices of block `blk`."""
    return np.arange(blk * block_size, (blk + 1) * block_size)


def _block_hit_rate(sel, active, block_size):
    """Fraction of active blocks with at least one selected variable."""
    if len(active) == 0:
        return 0.0
    sel_set = set(int(i) for i in sel)
    hits = sum(any(int(c) in sel_set for c in _block_cols(blk, block_size))
               for blk in active)
    return hits / len(active)


def _block_fdp(sel, active, G, block_size):
    """Block-level FDP = (null blocks hit) / max(total blocks hit, 1)."""
    sel_set = set(int(i) for i in sel)
    active_set = set(int(a) for a in active)
    total_hit = 0
    null_hit = 0
    for blk in range(G):
        if any(int(c) in sel_set for c in _block_cols(blk, block_size)):
            total_hit += 1
            if blk not in active_set:
                null_hit += 1
    return null_hit / max(total_hit, 1)


def _full_block_rate(sel, active, block_size):
    """Fraction of active blocks fully contained in the selection."""
    if len(active) == 0:
        return 0.0
    sel_set = set(int(i) for i in sel)
    fully = sum(all(int(c) in sel_set for c in _block_cols(blk, block_size))
                for blk in active)
    return fully / len(active)


def _block_purity_rate(groups_vec, G, block_size):
    """Fraction of the G true blocks whose members all share one groups_vec
    (discovered-cluster) label. Mirrors cpp compute_block_purity_rate; 1.0 for
    oracle runs by construction. `groups_vec` is the per-variable cluster label
    vector from the GVSSelectionResult (length p)."""
    gv = np.asarray(groups_vec)
    if gv.size == 0:
        return None
    pure = sum(bool(np.all(gv[k * block_size:(k + 1) * block_size] == gv[k * block_size]))
               for k in range(G))
    return pure / G


# ==============================================================================
# One method on one dataset
# ==============================================================================

def _run_block_method(dat, gvs_type, oracle, cv_seed, cfg):
    """Run one GVS method variant on a dataset and collect its metrics.

    `oracle` supplies the true block labels (prior_groups); otherwise groups are
    discovered from the data via HAC. `purity` and `M_found` are read from the
    GVSSelectionResult's `groups_vec` / `max_clusters` (available for both HAC and
    oracle runs).
    """
    en_solver = cfg["en_solver"] if gvs_type == "EN" else "TENET"
    groups = dat["prior_groups"] if oracle else None
    sel, res = build_gvs_selector(
        dat["X"], dat["y"], cfg["tFDR"], cfg["K"], gvs_type, cfg["corr_max"],
        hc_linkage="single", en_solver=en_solver, lambda2_method="CV_1SE_CCD",
        groups=groups, cv_seed=cv_seed, return_result=True)

    si = list(sel.selected_indices)
    ts = dat["true_support"].tolist()
    return dict(
        coord_tpr=compute_tpp(si, ts),
        coord_fdp=compute_fdp(si, ts),
        block_hit=_block_hit_rate(si, dat["active_blocks"], dat["block_size"]),
        block_fdp=_block_fdp(si, dat["active_blocks"], dat["G"], dat["block_size"]),
        full_block=_full_block_rate(si, dat["active_blocks"], dat["block_size"]),
        purity=_block_purity_rate(res.groups_vec, dat["G"], dat["block_size"]),
        T_stop=sel.T_stop,
        M_found=res.max_clusters,
    )


def _run_block_all_methods(dat, cv_seed, cfg):
    """Run all four method variants (M1-M4) on one dataset."""
    return dict(
        M1=_run_block_method(dat, "EN",  False, cv_seed, cfg),
        M2=_run_block_method(dat, "EN",  True,  cv_seed, cfg),
        M3=_run_block_method(dat, "IEN", False, cv_seed, cfg),
        M4=_run_block_method(dat, "IEN", True,  cv_seed, cfg),
    )


# ==============================================================================
# One (rho, snr) cell: num_MC trials, aggregated
# ==============================================================================

def _block_cell_worker(args):
    """One MC trial: generate data and run all four methods. Top-level so the
    spawn workers can import it from this module."""
    trial_seed, cfg, scenario, rho, snr = args
    dat = dgp_block_equicorr(
        n=cfg["n"], p=cfg["p"], G=cfg["G"], block_size=cfg["block_size"],
        n_active_blocks=cfg["n_active_blocks"], rho=rho, snr=snr,
        scenario=scenario, seed=trial_seed, b=cfg["b"])
    return _run_block_all_methods(dat, cv_seed=trial_seed, cfg=cfg)


def _agg_method(trials, mkey):
    """Average each metric across trials, treating None (NA) as missing."""
    out = {}
    for mt in _METRIC_NAMES:
        vals = np.array([np.nan if t[mkey][mt] is None else float(t[mkey][mt])
                         for t in trials])
        m = np.nanmean(vals) if np.any(~np.isnan(vals)) else np.nan
        out[mt] = m
    return out


def _run_block_cell(cfg, scenario, rho, snr, cell_base_seed, n_cores):
    tasks = [(cell_base_seed + mc, cfg, scenario, rho, snr)
             for mc in range(cfg["num_MC"])]
    if n_cores > 1:
        with multiprocessing.get_context("spawn").Pool(n_cores) as pool:
            trials = pool.map(_block_cell_worker, tasks)
    else:
        trials = [_block_cell_worker(t) for t in tasks]

    return {mk: _agg_method(trials, mk) for mk in ("M1", "M2", "M3", "M4")}


# ==============================================================================
# Grid table printer  (mirrors cpp print_block_grid_table)
# ==============================================================================

def _fmt(x, d):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return f"{'-':>9}"
    return f"{x:9.{d}f}"


def _print_block_grid(scenario, rho_grid, snr_grid, results, cfg):
    blk_metric_name = "full_blk" if scenario == "WHOLE" else "blk_hit"

    def get_blk(a):
        return a["full_block"] if scenario == "WHOLE" else a["block_hit"]

    methods = [("M1(EN-C)", "M1"), ("M2(EN-O)", "M2"),
               ("M3(IE-C)", "M3"), ("M4(IE-O)", "M4")]

    print("\n" + "=" * 90)
    print(f"  Block Bench MC: {scenario}  (MC={cfg['num_MC']}, tFDR={cfg['tFDR']:.2f},"
          f" n={cfg['n']}, p={cfg['p']}, G={cfg['G']}, blk={cfg['block_size']})")
    print("=" * 90)

    for ir, rho in enumerate(rho_grid):
        print(f"\n  rho = {rho:.2f}")
        print(f"  {'method':<12}{'SNR':>7}{'TPR':>9}{'FDR':>9}{'block_FDR':>9}"
              f"{blk_metric_name:>9}{'purity':>9}{'T_stop':>9}{'M_found':>9}")
        print("  " + "-" * (12 + 7 + 7 * 9))
        for isnr, snr in enumerate(snr_grid):
            cell = results[ir][isnr]
            for label, mk in methods:
                a = cell[mk]
                print(f"  {label:<12}{snr:7.2f}"
                      f"{_fmt(a['coord_tpr'], 4)}{_fmt(a['coord_fdp'], 4)}"
                      f"{_fmt(a['block_fdp'], 4)}{_fmt(get_blk(a), 4)}"
                      f"{_fmt(a['purity'], 4)}{_fmt(a['T_stop'], 1)}"
                      f"{_fmt(a['M_found'], 1)}")
            print()
    print("=" * 90 + "\n")


# ==============================================================================
# Benchmark driver for one truth scenario
# ==============================================================================

def _run_block_scenario(cfg, scenario):
    n_rho = len(cfg["rho_grid"])
    n_snr = len(cfg["snr_grid"])

    print("\n" + "=" * 90)
    print(f"Block Bench MC: {scenario}   n={cfg['n']} p={cfg['p']} G={cfg['G']}"
          f" blk={cfg['block_size']} active={cfg['n_active_blocks']} tFDR={cfg['tFDR']:.2f}"
          f" K={cfg['K']} MC={cfg['num_MC']} corr_max={cfg['corr_max']:.2f}")
    print("=" * 90)

    results = [[None] * n_snr for _ in range(n_rho)]
    for ir in range(n_rho):
        for isnr in range(n_snr):
            # Seeds staggered by grid-cell index so each cell uses a distinct band.
            cell_base_seed = cfg["base_seed"] + (ir * n_snr + isnr) * cfg["num_MC"]
            print(f"  [{scenario}] rho={cfg['rho_grid'][ir]:.2f}"
                  f" snr={cfg['snr_grid'][isnr]:.2f}  [{ir * n_snr + isnr + 1}/{n_rho * n_snr}]"
                  f"  {cfg['num_MC']} MC x 4 methods ...")
            results[ir][isnr] = _run_block_cell(
                cfg, scenario, cfg["rho_grid"][ir], cfg["snr_grid"][isnr],
                cell_base_seed, N_WORKERS)
    _print_block_grid(scenario, cfg["rho_grid"], cfg["snr_grid"], results, cfg)


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

def part_1_single_run():
    print("\n" + "=" * 70)
    print("Block-equicorrelated GVS  |  Part 1: Single-run")
    print(f"n={CFG['n']}, p={CFG['p']}, G={CFG['G']}, block_size={CFG['block_size']},"
          f" active blocks={CFG['n_active_blocks']}")
    print("Compares HAC-discovered vs oracle groups on one INDIVIDUAL-scenario draw.")
    print("=" * 70 + "\n")

    dat = dgp_block_equicorr(
        n=CFG["n"], p=CFG["p"], G=CFG["G"], block_size=CFG["block_size"],
        n_active_blocks=CFG["n_active_blocks"], rho=0.9, snr=2.0,
        scenario="INDIVIDUAL", seed=CFG["base_seed"], b=CFG["b"])
    print(f"[Part 1] Active blocks (IDs): {{{', '.join(map(str, dat['active_blocks']))}}}")
    print(f"[Part 1] True active variables: s = {dat['s']}"
          f"  |  sigma_y = {dat['sigma_y']:.4f}\n")

    labels = {"M1": "EN  + HAC ", "M2": "EN  + oracle",
              "M3": "IEN + HAC ", "M4": "IEN + oracle"}
    specs = {"M1": ("EN", False), "M2": ("EN", True),
             "M3": ("IEN", False), "M4": ("IEN", True)}
    print(f"  {'method':<14}{'TPR':>8}{'FDR':>8}{'block_FDR':>10}{'blk_hit':>9}"
          f"{'purity':>9}{'T_stop':>9}{'M_found':>9}")
    print("  " + "-" * 76)
    for mk in ("M1", "M2", "M3", "M4"):
        gvs_type, oracle = specs[mk]
        r = _run_block_method(dat, gvs_type, oracle, cv_seed=CFG["base_seed"], cfg=CFG)
        print(f"  {labels[mk]:<14}{r['coord_tpr']:8.4f}{r['coord_fdp']:8.4f}"
              f"{r['block_fdp']:10.4f}{r['block_hit']:9.4f}{_fmt(r['purity'], 4)}"
              f"{r['T_stop']:9.1f}{_fmt(r['M_found'], 1)}")
    print()


# ==============================================================================
# Part 2: Full benchmark — all three truth scenarios
# ==============================================================================

def part_2_full_benchmark():
    print("\n" + "=" * 90)
    print("Block-equicorrelated benchmark")
    print("methods = {EN clustered, EN oracle, IEN clustered, IEN oracle}")
    print("scenarios = {INDIVIDUAL, REPRESENTATIVE, WHOLE}")
    print("=" * 90)

    for scenario in ("INDIVIDUAL", "REPRESENTATIVE", "WHOLE"):
        _run_block_scenario(CFG, scenario)


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    if RUN_PART_1:
        part_1_single_run()
    if RUN_PART_2:
        part_2_full_benchmark()
    print("\nBlock benchmark GVS demo complete.")
