# ==============================================================================
# demo_trex_gvs_07.py
# ==============================================================================
#
# T-Rex+GVS on the ARMA mixed-structure blocks DGP.
#
# Four blocks (sizes {20, 50, 80, 65}) each following a different time-series
# model:
#   Block 0 = AR(1)      [size 20]  — active
#   Block 1 = MA(3)      [size 50]  — active
#   Block 2 = ARMA(2,1)  [size 80]  — active   (s = 150)
#   Block 3 = AR(1) trap [size 65]  — inactive
# Shuffled into p = 500 columns with white-noise gaps. The AR component is
# governed by ar_coef; ma_coefs = (0.5, 0.3, 0.1) is fixed. This is the most
# heterogeneous within-block structure in the suite.
#
# Methods compared: EN (TENET), EN+AUG (TENET_AUG), IEN.
# The swept correlation knob here is ar_coef (not rho).
#
# Parts:
#   Part 1: Single-run demo — MA(3) lag-structure check + block layout + one run.
#   Part 2: SNR sweep       — fixed ar_coef = 0.8.
#   Part 3: ar_coef sweep   — fixed SNR = 2.0.
#   Part 4: 2-D SNR x ar_coef grid.
#
# Mirrors R/trex_selector_methods/trex_gvs/demo_trex_gvs_07. Part 1 runs by
# default; the heavy MC parts are opt-in via RUN_PART_*.
# Indices are 0-based (Python convention); the R counterpart is 1-based.
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

from trex_gvs_dgps import dgp_arma_blocks_snr
from gvs_sim_common import (N_CORES, build_gvs_selector, run_mc_gvs,
                            print_table, print_matrix, print_single_run_result)

N_WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else N_CORES

CFG = dict(
    n=200, p=500,
    ar_coef=0.8,            # fixed AR parameter for the SNR sweep
    snr=2.0,                # fixed SNR for the ar_coef sweep
    K=20, tFDR=0.1,
    seed=2026, num_MC=200,
    corr_max=0.98, hc_dist="single",
)

GVS_METHODS = [
    dict(label="EN",     gvs_type="EN",  en_solver="TENET"),
    dict(label="EN+AUG", gvs_type="EN",  en_solver="TENET_AUG"),
    dict(label="IEN",    gvs_type="IEN", en_solver="TENET"),
]

RUN_PART_1 = True
RUN_PART_2 = False
RUN_PART_3 = False
RUN_PART_4 = False


def _arma_point(method, snr, ar_coef):
    """Run one MC point on the ARMA-blocks DGP for (method, snr, ar_coef)."""
    return run_mc_gvs(
        "dgp_arma_blocks_snr",
        dict(n=CFG["n"], p=CFG["p"], snr=snr, ar_coef=ar_coef),
        CFG["tFDR"], CFG["K"], CFG["num_MC"], CFG["seed"],
        method["gvs_type"], CFG["corr_max"], hc_linkage=CFG["hc_dist"],
        en_solver=method["en_solver"], n_cores=N_WORKERS)


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

def part_1_single_run():
    print("\n" + "=" * 70)
    print("ARMA-Blocks GVS  |  Part 1: Single-run")
    print(f"n={CFG['n']},  p={CFG['p']},  s=150  (blocks 20+50+80 active; 65 trap)")
    print(f"ar_coef={CFG['ar_coef']:.2f},  SNR={CFG['snr']:.2f},  tFDR={CFG['tFDR']:.2f}")
    print("Block types: AR(1) | MA(3) | ARMA(2,1) | AR(1) trap")
    print("=" * 70 + "\n")

    print("[Part 1] Generating ARMA-blocks data ...")
    dat = dgp_arma_blocks_snr(n=CFG["n"], p=CFG["p"], snr=CFG["snr"],
                              ar_coef=CFG["ar_coef"], seed=CFG["seed"])
    print(f"[Part 1] Active: {dat['s']}  |  sigma_y = {dat['sigma_y']:.4f}\n")

    print("[Part 1] Block layout (shuffled order):")
    block_sizes = [20, 50, 80, 65]
    block_types = ["AR(1)", "MA(3)", "ARMA(2,1)", "AR(1) trap"]
    active_label = ["active", "active", "active", "trap"]
    for b in range(4):
        idx = dat["block_indices"][b]
        print(f"  Block {b} [size={block_sizes[b]:2d}, {block_types[b]}, {active_label[b]}]:"
              f" columns {int(idx.min())} - {int(idx.max())}")
    print(f"  Block order (ID sequence placed): {{{', '.join(map(str, dat['block_order']))}}}\n")

    # MA(3) lag structure: positive short-lag correlation, ~zero beyond lag 3.
    ma3 = dat["block_indices"][1]
    X = dat["X"]
    lag1 = np.corrcoef(X[:, ma3[0]], X[:, ma3[1]])[0, 1]
    lag4 = np.corrcoef(X[:, ma3[0]], X[:, ma3[4]])[0, 1]
    print(f"[Part 1] MA(3) Block 1 - Lag 1 correlation (expect > 0): {lag1:.4f}")
    print(f"[Part 1] MA(3) Block 1 - Lag 4 (expect ~0):              {lag4:.4f}\n")

    for m in GVS_METHODS:
        print(f"[Part 1] Running trex+GVS ({m['label']}) ...\n")
        sel = build_gvs_selector(dat["X"], dat["y"], CFG["tFDR"], CFG["K"],
                                 m["gvs_type"], CFG["corr_max"], CFG["hc_dist"],
                                 m["en_solver"])
        print_single_run_result(f"Part 1 - ARMA-Blocks GVS  [{m['label']}]",
                                 dat, sel, CFG["tFDR"])


# ==============================================================================
# Part 2: MC simulation - SNR sweep  (fixed ar_coef = 0.8)
# ==============================================================================

def part_2_snr_sweep():
    snr_grid = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

    print("=" * 70)
    print("ARMA-Blocks GVS MC - Part 2: SNR sweep")
    print(f"tFDR={CFG['tFDR']:.2f}  ar_coef={CFG['ar_coef']:.2f}  n={CFG['n']}"
          f"  p={CFG['p']}  {CFG['num_MC']} MC")
    print("=" * 70 + "\n")

    for m in GVS_METHODS:
        print(f"\n  [{m['label']}] SNR sweep, {CFG['num_MC']} MC per point ...")
        rows = []
        for s in snr_grid:
            print(f"    SNR = {s:.2f} ...")
            rows.append({"SNR": s, **_arma_point(m, s, CFG["ar_coef"])})
        print(f"\n  Results - {m['label']}:")
        print_table(rows, "SNR")


# ==============================================================================
# Part 3: MC simulation - ar_coef sweep  (fixed SNR = 2.0)
# ==============================================================================

def part_3_ar_coef_sweep():
    ar_coef_grid = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95]

    print("=" * 70)
    print("ARMA-Blocks GVS MC - Part 3: ar_coef sweep")
    print(f"tFDR={CFG['tFDR']:.2f}  SNR={CFG['snr']:.2f}  n={CFG['n']}"
          f"  p={CFG['p']}  {CFG['num_MC']} MC")
    print("=" * 70 + "\n")

    for m in GVS_METHODS:
        print(f"\n  [{m['label']}] ar_coef sweep, {CFG['num_MC']} MC per point ...")
        rows = []
        for ar in ar_coef_grid:
            print(f"    ar_coef = {ar:.2f} ...")
            rows.append({"ar_coef": ar, **_arma_point(m, CFG["snr"], ar)})
        print(f"\n  Results - {m['label']}:")
        print_table(rows, "ar_coef", value_fmt="{:<8.2f}", width=8)


# ==============================================================================
# Part 4: 2-D phase-transition sweep - SNR x ar_coef
# ==============================================================================

def part_4_snr_ar_grid():
    snr_grid = [0.2, 0.5, 1.0, 2.0, 5.0]
    ar_coef_grid = [0.30, 0.50, 0.70, 0.80, 0.90, 0.95]
    row_labels = [f"snr={s:.2f}" for s in snr_grid]
    col_labels = [f"ar={a:.2f}" for a in ar_coef_grid]

    print("=" * 70)
    print("ARMA-Blocks GVS MC  |  Part 4: SNR x ar_coef grid")
    print(f"n={CFG['n']}  p={CFG['p']}  {CFG['num_MC']} MC  tFDR={CFG['tFDR']:.2f}"
          f"  corr_max={CFG['corr_max']:.2f}")
    print("=" * 70 + "\n")

    for m in GVS_METHODS:
        mat_tpp = np.full((len(snr_grid), len(ar_coef_grid)), np.nan)
        mat_fdp = np.full_like(mat_tpp, np.nan)
        for i, s in enumerate(snr_grid):
            for j, ar in enumerate(ar_coef_grid):
                print(f"  [{m['label']}] SNR={s:.2f} ar_coef={ar:.2f} ...")
                r = _arma_point(m, s, ar)
                mat_tpp[i, j] = r["mean_TPP"]
                mat_fdp[i, j] = r["mean_FDP"]
        print_matrix(mat_tpp, f"mean_TPP [{m['label']}] (rows: SNR, cols: ar_coef)",
                     row_labels, col_labels)
        print_matrix(mat_fdp, f"mean_FDP [{m['label']}] (rows: SNR, cols: ar_coef)",
                     row_labels, col_labels)


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    if RUN_PART_1:
        part_1_single_run()
    if RUN_PART_2:
        part_2_snr_sweep()
    if RUN_PART_3:
        part_3_ar_coef_sweep()
    if RUN_PART_4:
        part_4_snr_ar_grid()
    print("\nARMA-blocks GVS demo complete.")
