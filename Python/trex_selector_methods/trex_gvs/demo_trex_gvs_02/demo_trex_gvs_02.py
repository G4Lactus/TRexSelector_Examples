# ==============================================================================
# demo_trex_gvs_02.py
# ==============================================================================
#
# T-Rex+GVS on the scattered-grouped DGP.
#
# Three equicorrelated groups of 50 variables (all active, s = 150) whose members
# are randomly scattered across the p = 500 columns instead of being contiguous.
# This isolates the effect of column layout from correlation strength: the groups
# are identical to demo 01, only their positions are permuted.
#
# Methods compared: EN (TENET), EN+AUG (TENET_AUG), IEN.
# Within-group correlation is controlled through sd_x, with rho = 1/(1 + sd_x^2).
#
# Parts:
#   Part 1: Single-run demo  — one EN/EN+AUG/IEN run.
#   Part 2: SNR sweep        — fixed sd_x = sqrt(0.01) (rho ~ 0.99).
#   Part 3: rho sweep        — fixed SNR = 2.0 (sd_x derived from rho).
#   Part 4: 2-D SNR x rho grid.
#
# Mirrors R/trex_selector_methods/trex_gvs/demo_trex_gvs_02. Part 1 runs by
# default; the heavy MC parts (200 MC per grid point) are opt-in via RUN_PART_*.
# Indices are 0-based (Python convention); the R counterpart is 1-based.
# ==============================================================================

import sys
import os
import math

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

from trex_gvs_dgps import dgp_scattered_grouped_snr
from gvs_sim_common import (N_CORES, build_gvs_selector, run_mc_gvs,
                            print_table, print_matrix, print_single_run_result)

N_WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else N_CORES

CFG = dict(
    n=200, p=500,
    n_groups=3, group_size=50,
    sd_x=math.sqrt(0.01),   # rho ~ 0.99
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


def _scattered_point(method, snr, sd_x):
    """Run one MC point on the scattered-grouped DGP for (method, snr, sd_x)."""
    return run_mc_gvs(
        "dgp_scattered_grouped_snr",
        dict(n=CFG["n"], p=CFG["p"], snr=snr, sd_x=sd_x,
             n_groups=CFG["n_groups"], group_size=CFG["group_size"]),
        CFG["tFDR"], CFG["K"], CFG["num_MC"], CFG["seed"],
        method["gvs_type"], CFG["corr_max"], hc_linkage=CFG["hc_dist"],
        en_solver=method["en_solver"], n_cores=N_WORKERS)


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

def part_1_single_run():
    print("\n" + "=" * 70)
    print("Scattered-Grouped GVS  |  Part 1: Single-run")
    print(f"n={CFG['n']},  p={CFG['p']},  n_groups={CFG['n_groups']},"
          f"  group_size={CFG['group_size']}")
    print(f"SNR={CFG['snr']:.2f},  sd_x={CFG['sd_x']:.4f},  tFDR={CFG['tFDR']:.2f}")
    print("=" * 70 + "\n")

    print("[Part 1] Generating scattered-grouped data ...")
    dat = dgp_scattered_grouped_snr(
        n=CFG["n"], p=CFG["p"], snr=CFG["snr"], sd_x=CFG["sd_x"],
        n_groups=CFG["n_groups"], group_size=CFG["group_size"], seed=CFG["seed"])
    print(f"[Part 1] Active variables: {dat['s']}  |  sigma_y = {dat['sigma_y']:.4f}\n")

    for m in GVS_METHODS:
        print(f"[Part 1] Running trex+GVS ({m['label']}) ...\n")
        sel = build_gvs_selector(dat["X"], dat["y"], CFG["tFDR"], CFG["K"],
                                 m["gvs_type"], CFG["corr_max"], CFG["hc_dist"],
                                 m["en_solver"])
        print_single_run_result(f"Part 1 - Scattered-Grouped GVS  [{m['label']}]",
                                 dat, sel, CFG["tFDR"])


# ==============================================================================
# Part 2: MC simulation - SNR sweep  (fixed sd_x = sqrt(0.01), rho ~ 0.99)
# ==============================================================================

def part_2_snr_sweep():
    snr_grid = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

    print("=" * 70)
    print("Scattered-Grouped GVS MC - Part 2: SNR sweep")
    print(f"tFDR={CFG['tFDR']:.2f}  sd_x={CFG['sd_x']:.4f}"
          f" (rho~{1 / (1 + CFG['sd_x'] ** 2):.2f})  n={CFG['n']}  p={CFG['p']}"
          f"  {CFG['num_MC']} MC")
    print("=" * 70 + "\n")

    for m in GVS_METHODS:
        print(f"\n  [{m['label']}] SNR sweep, {CFG['num_MC']} MC per point ...")
        rows = []
        for s in snr_grid:
            print(f"    SNR = {s:.2f} ...")
            rows.append({"SNR": s, **_scattered_point(m, s, CFG["sd_x"])})
        print(f"\n  Results - {m['label']}:")
        print_table(rows, "SNR")


# ==============================================================================
# Part 3: MC simulation - rho sweep  (fixed SNR = 2.0)
# ==============================================================================

def part_3_rho_sweep():
    rho_grid = [0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99]

    print("=" * 70)
    print("Scattered-Grouped GVS MC - Part 3: rho sweep")
    print(f"tFDR={CFG['tFDR']:.2f}  SNR={CFG['snr']:.2f}  n={CFG['n']}"
          f"  p={CFG['p']}  {CFG['num_MC']} MC")
    print("=" * 70 + "\n")

    for m in GVS_METHODS:
        print(f"\n  [{m['label']}] rho sweep, {CFG['num_MC']} MC per point ...")
        rows = []
        for rho in rho_grid:
            sd_x = math.sqrt((1 - rho) / rho)
            print(f"    rho = {rho:.2f}  (sd_x = {sd_x:.4f}) ...")
            rows.append({"rho": rho, **_scattered_point(m, CFG["snr"], sd_x)})
        print(f"\n  Results - {m['label']}:")
        print_table(rows, "rho", value_fmt="{:<6.2f}", width=6)


# ==============================================================================
# Part 4: 2-D phase-transition sweep - SNR x rho
# ==============================================================================

def part_4_snr_rho_grid():
    snr_grid = [0.2, 0.5, 1.0, 2.0, 5.0]
    rho_grid = [0.30, 0.50, 0.70, 0.90, 0.95, 0.99]
    row_labels = [f"snr={s:.2f}" for s in snr_grid]
    col_labels = [f"rho={r:.2f}" for r in rho_grid]

    print("=" * 70)
    print("Scattered-Grouped GVS MC  |  Part 4: SNR x rho grid")
    print(f"n={CFG['n']}  p={CFG['p']}  {CFG['num_MC']} MC  tFDR={CFG['tFDR']:.2f}"
          f"  corr_max={CFG['corr_max']:.2f}")
    print("=" * 70 + "\n")

    for m in GVS_METHODS:
        mat_tpp = np.full((len(snr_grid), len(rho_grid)), np.nan)
        mat_fdp = np.full_like(mat_tpp, np.nan)
        for i, s in enumerate(snr_grid):
            for j, rho in enumerate(rho_grid):
                sd_x = math.sqrt((1 - rho) / rho)
                print(f"  [{m['label']}] SNR={s:.2f} rho={rho:.2f} ...")
                r = _scattered_point(m, s, sd_x)
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
    print("\nScattered-grouped GVS demo complete.")
