# ==============================================================================
# demo_trex_gvs_06.py
# ==============================================================================
#
# T-Rex+GVS on the block-structured AR(1) DGP.
#
# Same 4-block geometry as demo 03 (sizes {20, 50, 80, 65}: 3 active, s = 150,
# plus 1 inactive trap, shuffled into p = 500 with white-noise gaps), but the
# within-block correlation is AR(1)/Toeplitz: Cor(X_j, X_k) = rho^|j-k|, so
# correlation decays with column distance instead of being equicorrelated.
# Robustness to non-equicorrelated within-block structure.
#
# Methods compared: EN (TENET), EN+AUG (TENET_AUG), IEN.
# Within-block AR(1) parameter rho is passed directly to the DGP.
#
# Parts:
#   Part 1: Single-run demo — AR(1) decay check + block layout + one EN/EN+AUG/IEN run.
#   Part 2: SNR sweep       — fixed rho = 0.85.
#   Part 3: rho sweep       — fixed SNR = 2.0.
#   Part 4: 2-D SNR x rho grid.
#
# Mirrors R/trex_selector_methods/trex_gvs/demo_trex_gvs_06. Part 1 runs by
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

from trex_gvs_dgps import dgp_ar1_blocks_snr
from gvs_sim_common import (N_CORES, build_gvs_selector, run_mc_gvs,
                            print_table, print_matrix, print_single_run_result)

N_WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else N_CORES

CFG = dict(
    n=200, p=500,
    rho=0.85,               # fixed within-block AR(1) parameter for the SNR sweep
    snr=2.0,                # fixed SNR for the rho sweep
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


def _ar1_point(method, snr, rho):
    """Run one MC point on the AR(1)-blocks DGP for (method, snr, rho)."""
    return run_mc_gvs(
        "dgp_ar1_blocks_snr", dict(n=CFG["n"], p=CFG["p"], snr=snr, rho=rho),
        CFG["tFDR"], CFG["K"], CFG["num_MC"], CFG["seed"],
        method["gvs_type"], CFG["corr_max"], hc_linkage=CFG["hc_dist"],
        en_solver=method["en_solver"], n_cores=N_WORKERS)


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

def part_1_single_run():
    print("\n" + "=" * 70)
    print("AR(1)-Blocks GVS  |  Part 1: Single-run")
    print(f"n={CFG['n']},  p={CFG['p']},  s=150  (blocks 20+50+80 active; 65 trap)")
    print(f"rho={CFG['rho']:.2f},  SNR={CFG['snr']:.2f},  tFDR={CFG['tFDR']:.2f}")
    print("=" * 70 + "\n")

    print("[Part 1] Generating AR(1)-blocks data ...")
    dat = dgp_ar1_blocks_snr(n=CFG["n"], p=CFG["p"], snr=CFG["snr"],
                             rho=CFG["rho"], seed=CFG["seed"])
    print(f"[Part 1] Active: {dat['s']}  |  sigma_y = {dat['sigma_y']:.4f}\n")

    print("[Part 1] Block layout (shuffled order):")
    block_sizes = [20, 50, 80, 65]
    active_label = ["active", "active", "active", "trap"]
    for b in range(4):
        idx = dat["block_indices"][b]
        print(f"  Block {b} [size={block_sizes[b]:2d}, {active_label[b]}]:"
              f" columns {int(idx.min())} - {int(idx.max())}")
    print(f"  Block order (ID sequence placed): {{{', '.join(map(str, dat['block_order']))}}}\n")

    # AR(1) decay: empirical vs theoretical rho^k within block 0.
    b1 = dat["block_indices"][0]
    X = dat["X"]
    rho1 = np.corrcoef(X[:, b1[0]], X[:, b1[1]])[0, 1]
    rho2 = np.corrcoef(X[:, b1[0]], X[:, b1[2]])[0, 1]
    rho3 = np.corrcoef(X[:, b1[0]], X[:, b1[3]])[0, 1]
    print("[Part 1] AR(1) decay (theoretical rho^k vs empirical):")
    print(f"  Lag 1: {rho1:.4f} (target rho   = {CFG['rho']:.4f})")
    print(f"  Lag 2: {rho2:.4f} (target rho^2 = {CFG['rho'] ** 2:.4f})")
    print(f"  Lag 3: {rho3:.4f} (target rho^3 = {CFG['rho'] ** 3:.4f})\n")

    for m in GVS_METHODS:
        print(f"[Part 1] Running trex+GVS ({m['label']}) ...\n")
        sel = build_gvs_selector(dat["X"], dat["y"], CFG["tFDR"], CFG["K"],
                                 m["gvs_type"], CFG["corr_max"], CFG["hc_dist"],
                                 m["en_solver"])
        print_single_run_result(f"Part 1 - AR(1)-Blocks GVS  [{m['label']}]",
                                 dat, sel, CFG["tFDR"])


# ==============================================================================
# Part 2: MC simulation - SNR sweep  (fixed rho = 0.85)
# ==============================================================================

def part_2_snr_sweep():
    snr_grid = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

    print("=" * 70)
    print("AR(1)-Blocks GVS MC - Part 2: SNR sweep")
    print(f"tFDR={CFG['tFDR']:.2f}  rho={CFG['rho']:.2f}  n={CFG['n']}"
          f"  p={CFG['p']}  {CFG['num_MC']} MC")
    print("=" * 70 + "\n")

    for m in GVS_METHODS:
        print(f"\n  [{m['label']}] SNR sweep, {CFG['num_MC']} MC per point ...")
        rows = []
        for s in snr_grid:
            print(f"    SNR = {s:.2f} ...")
            rows.append({"SNR": s, **_ar1_point(m, s, CFG["rho"])})
        print(f"\n  Results - {m['label']}:")
        print_table(rows, "SNR")


# ==============================================================================
# Part 3: MC simulation - rho sweep  (fixed SNR = 2.0)
# ==============================================================================

def part_3_rho_sweep():
    rho_grid = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]

    print("=" * 70)
    print("AR(1)-Blocks GVS MC - Part 3: rho sweep")
    print(f"tFDR={CFG['tFDR']:.2f}  SNR={CFG['snr']:.2f}  n={CFG['n']}"
          f"  p={CFG['p']}  {CFG['num_MC']} MC")
    print("=" * 70 + "\n")

    for m in GVS_METHODS:
        print(f"\n  [{m['label']}] rho sweep, {CFG['num_MC']} MC per point ...")
        rows = []
        for rho in rho_grid:
            print(f"    rho = {rho:.2f} ...")
            rows.append({"rho": rho, **_ar1_point(m, CFG["snr"], rho)})
        print(f"\n  Results - {m['label']}:")
        print_table(rows, "rho", value_fmt="{:<6.2f}", width=6)


# ==============================================================================
# Part 4: 2-D phase-transition sweep - SNR x rho
# ==============================================================================

def part_4_snr_rho_grid():
    snr_grid = [0.2, 0.5, 1.0, 2.0, 5.0]
    rho_grid = [0.30, 0.50, 0.70, 0.85, 0.90, 0.99]
    row_labels = [f"snr={s:.2f}" for s in snr_grid]
    col_labels = [f"rho={r:.2f}" for r in rho_grid]

    print("=" * 70)
    print("AR(1)-Blocks GVS MC  |  Part 4: SNR x rho grid")
    print(f"n={CFG['n']}  p={CFG['p']}  {CFG['num_MC']} MC  tFDR={CFG['tFDR']:.2f}"
          f"  corr_max={CFG['corr_max']:.2f}")
    print("=" * 70 + "\n")

    for m in GVS_METHODS:
        mat_tpp = np.full((len(snr_grid), len(rho_grid)), np.nan)
        mat_fdp = np.full_like(mat_tpp, np.nan)
        for i, s in enumerate(snr_grid):
            for j, rho in enumerate(rho_grid):
                print(f"  [{m['label']}] SNR={s:.2f} rho={rho:.2f} ...")
                r = _ar1_point(m, s, rho)
                mat_tpp[i, j] = r["mean_TPP"]
                mat_fdp[i, j] = r["mean_FDP"]
        print_matrix(mat_tpp, f"mean_TPP [{m['label']}] (rows: SNR, cols: rho)",
                     row_labels, col_labels)
        print_matrix(mat_fdp, f"mean_FDP [{m['label']}] (rows: SNR, cols: rho)",
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
        part_3_rho_sweep()
    if RUN_PART_4:
        part_4_snr_rho_grid()
    print("\nAR(1)-blocks GVS demo complete.")
