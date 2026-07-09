# ==============================================================================
# demo_trex_gvs_03.py
# ==============================================================================
#
# T-Rex+GVS on the mixed-blocks DGP.
#
# Four unequal blocks of sizes {20, 50, 80, 65}: three active (beta = 3, s = 150)
# plus one inactive equicorrelated trap block, placed in a random order separated
# by white-noise gaps. A more realistic grouped design than the contiguous Hastie
# blocks of demo 01.
#
# Methods compared: EN (TENET), EN+AUG (TENET_AUG), IEN.
#
# Following cpp/R demo 03, two presets are swept (each = SNR + rho + 2-D grid):
#   * Mixed-Blocks        — default within-block noise sd_x = sqrt(0.01) (rho ~ 0.99).
#   * Mixed-Blocks-Rho075 — fixed rho = 0.75 for the SNR sweep; a wider rho grid.
#
# Parts:
#   Part 1: Single-run demo — block-layout diagnostics + one EN/EN+AUG/IEN run.
#   Part 2: Mixed-Blocks preset        — SNR sweep + rho sweep + 2-D grid.
#   Part 3: Mixed-Blocks-Rho075 preset — SNR sweep + rho sweep + 2-D grid.
#
# Mirrors R/trex_selector_methods/trex_gvs/demo_trex_gvs_03. Part 1 runs by
# default; the heavy MC presets are opt-in via RUN_PART_*.
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

from trex_gvs_dgps import dgp_mixed_blocks_snr
from gvs_sim_common import (N_CORES, build_gvs_selector, run_mc_gvs,
                            print_table, print_matrix, print_single_run_result)

N_WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else N_CORES

CFG = dict(
    n=200, p=500,
    sd_x=math.sqrt(0.01),   # default within-block noise sd (rho ~ 0.99)
    snr=2.0,                # fixed SNR for the rho sweeps
    K=20, tFDR=0.1,
    seed=2026, num_MC=200,
    corr_max=0.98, hc_dist="single",
)

GVS_METHODS = [
    dict(label="EN",     gvs_type="EN",  en_solver="TENET"),
    dict(label="EN+AUG", gvs_type="EN",  en_solver="TENET_AUG"),
    dict(label="IEN",    gvs_type="IEN", en_solver="TENET"),
]

# Two presets swept in Parts 2 and 3 (mirrors cpp/R demo 03).
PRESETS = dict(
    default_sd_x=dict(
        tag="Mixed-Blocks",
        snr_fixed_sdx=math.sqrt(0.01),                    # rho ~ 0.99
        snr_note="sd_x=sqrt(0.01)",
        snr_grid=[0.1, 0.2, 0.5, 1.0, 2.0, 5.0],
        rho_grid=[0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99],
        snr_grid_2d=[0.2, 0.5, 1.0, 2.0, 5.0],
        rho_grid_2d=[0.30, 0.50, 0.70, 0.90, 0.95, 0.99],
    ),
    fixed_rho_075=dict(
        tag="Mixed-Blocks-Rho075",
        snr_fixed_sdx=math.sqrt((1 - 0.75) / 0.75),       # rho = 0.75
        snr_note="rho=0.75",
        snr_grid=[0.1, 0.2, 0.5, 1.0, 2.0, 5.0],
        rho_grid=[0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99],
        snr_grid_2d=[0.2, 0.5, 1.0, 2.0, 5.0],
        rho_grid_2d=[0.30, 0.50, 0.70, 0.85, 0.90, 0.99],
    ),
)

RUN_PART_1 = True
RUN_PART_2 = False
RUN_PART_3 = False


def _mixed_point(method, snr, sd_x):
    """Run one MC point on the mixed-blocks DGP for (method, snr, sd_x)."""
    return run_mc_gvs(
        "dgp_mixed_blocks_snr", dict(n=CFG["n"], p=CFG["p"], snr=snr, sd_x=sd_x),
        CFG["tFDR"], CFG["K"], CFG["num_MC"], CFG["seed"],
        method["gvs_type"], CFG["corr_max"], hc_linkage=CFG["hc_dist"],
        en_solver=method["en_solver"], n_cores=N_WORKERS)


def _run_mixed_benchmark(preset):
    """Run the SNR + rho + 2-D benchmark for one preset."""
    print("=" * 70)
    print(f"Mixed-Blocks GVS MC  |  Preset: {preset['tag']}")
    print(f"n={CFG['n']}  p={CFG['p']}  {CFG['num_MC']} MC  tFDR={CFG['tFDR']:.2f}"
          f"  corr_max={CFG['corr_max']:.2f}")
    print("=" * 70 + "\n")

    # --- SNR sweep (fixed sd_x) ---
    print(f"  [{preset['tag']}] SNR sweep  (fixed {preset['snr_note']})")
    for m in GVS_METHODS:
        rows = []
        for s in preset["snr_grid"]:
            print(f"    [{m['label']}] SNR = {s:.2f} ...")
            rows.append({"SNR": s, **_mixed_point(m, s, preset["snr_fixed_sdx"])})
        print(f"\n  SNR sweep - {m['label']}:")
        print_table(rows, "SNR")

    # --- rho sweep (fixed SNR) ---
    print(f"  [{preset['tag']}] rho sweep  (fixed SNR = {CFG['snr']:.2f})")
    for m in GVS_METHODS:
        rows = []
        for rho in preset["rho_grid"]:
            sd_x = math.sqrt((1 - rho) / rho)
            print(f"    [{m['label']}] rho = {rho:.2f}  (sd_x = {sd_x:.4f}) ...")
            rows.append({"rho": rho, **_mixed_point(m, CFG["snr"], sd_x)})
        print(f"\n  rho sweep - {m['label']}:")
        print_table(rows, "rho", value_fmt="{:<6.2f}", width=6)

    # --- 2-D SNR x rho grid ---
    print(f"  [{preset['tag']}] 2-D SNR x rho grid")
    row_labels = [f"snr={s:.2f}" for s in preset["snr_grid_2d"]]
    col_labels = [f"rho={r:.2f}" for r in preset["rho_grid_2d"]]
    for m in GVS_METHODS:
        mat_tpp = np.full((len(preset["snr_grid_2d"]), len(preset["rho_grid_2d"])), np.nan)
        mat_fdp = np.full_like(mat_tpp, np.nan)
        for i, s in enumerate(preset["snr_grid_2d"]):
            for j, rho in enumerate(preset["rho_grid_2d"]):
                sd_x = math.sqrt((1 - rho) / rho)
                print(f"    [{m['label']}] SNR={s:.2f} rho={rho:.2f} ...")
                r = _mixed_point(m, s, sd_x)
                mat_tpp[i, j] = r["mean_TPP"]
                mat_fdp[i, j] = r["mean_FDP"]
        print_matrix(mat_tpp, f"mean_TPP [{m['label']}] (rows: SNR, cols: rho)",
                     row_labels, col_labels)
        print_matrix(mat_fdp, f"mean_FDP [{m['label']}] (rows: SNR, cols: rho)",
                     row_labels, col_labels)


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

def part_1_single_run():
    print("\n" + "=" * 70)
    print("Mixed-Blocks GVS  |  Part 1: Single-run")
    print(f"n={CFG['n']},  p={CFG['p']},  active blocks={{20,50,80}}, inactive={{65}}")
    print(f"SNR={CFG['snr']:.2f},  sd_x={CFG['sd_x']:.4f},  tFDR={CFG['tFDR']:.2f}")
    print("=" * 70 + "\n")

    print("[Part 1] Generating mixed-blocks data ...")
    dat = dgp_mixed_blocks_snr(n=CFG["n"], p=CFG["p"], snr=CFG["snr"],
                               sd_x=CFG["sd_x"], seed=CFG["seed"])
    print(f"[Part 1] Active variables: {dat['s']}  |  sigma_y = {dat['sigma_y']:.4f}")
    print(f"[Part 1] Block order (IDs): {{{', '.join(map(str, dat['block_order']))}}}")
    for b in range(4):
        idx = dat["block_indices"][b]
        lbl = "ACTIVE  " if dat["is_active"][b] else "inactive"
        print(f"[Part 1]   Block {b} (size {int(dat['block_sizes'][b]):2d}, {lbl}):"
              f" cols {int(idx.min()):3d} - {int(idx.max()):3d}")
    print()

    for m in GVS_METHODS:
        print(f"[Part 1] Running trex+GVS ({m['label']}) ...\n")
        sel = build_gvs_selector(dat["X"], dat["y"], CFG["tFDR"], CFG["K"],
                                 m["gvs_type"], CFG["corr_max"], CFG["hc_dist"],
                                 m["en_solver"])
        print_single_run_result(f"Part 1 - Mixed-Blocks GVS  [{m['label']}]",
                                 dat, sel, CFG["tFDR"])


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    if RUN_PART_1:
        part_1_single_run()
    if RUN_PART_2:
        _run_mixed_benchmark(PRESETS["default_sd_x"])
    if RUN_PART_3:
        _run_mixed_benchmark(PRESETS["fixed_rho_075"])
    print("\nMixed-blocks GVS demo complete.")
