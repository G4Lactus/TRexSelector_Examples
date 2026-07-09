# ==============================================================================
# demo_trex_da_01_ar1.py
# ==============================================================================
#
# DA-TRex demo and MC simulation for the AR(1) Toeplitz DGP.
#
#   Part 1: single-run demo.
#   Part 2: MC SNR sweep — Random vs CappedSpread support.
#   Part 3: MC rho sweep — Random vs CappedSpread support.
#   Part 4: 2D Gap x Rho sweep — probes the kappa boundary where gap < kappa
#           collapses TPP.
#
# Mirrors R/trex_selector_methods/trex_da/demo_trex_da_01_ar1.R and
# cpp/.../demo_trex_da_01_mc_sim_ar1. The heavy MC parts are opt-in via the
# RUN_PART_* flags (Part 1 runs by default).
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

from dgp_generators import dgp_ar1_snr
from da_sim_common import (make_support_capped_spread, make_support_random,
                           build_da_selector, run_mc, print_table, print_matrix,
                           compute_fdp, compute_tpp)

RUN_PART_1 = True
RUN_PART_2 = False
RUN_PART_3 = False
RUN_PART_4 = False

# Part 1 parameters
PARAMS = dict(n=150, p=500, s=5, amplitude=3.0, snr=2.0, K=20, tFDR=0.2,
              seed=2026, rho=0.7, max_gap=15)

# Parts 2-4 Monte Carlo base parameters
MC = dict(n=300, p=1000, s=10, amplitude=3.0, rho=0.7, snr=2.0, tFDR=0.2,
          K=20, num_MC=200, seed=2026, max_gap=20)

_DA_EXTRA = {"rho_thr_DA": 0.02}   # trex+DA+AR1: auto cor_coef, rho threshold 0.02


# ==============================================================================
# Part 1: single-run demo
# ==============================================================================

def part_1_single_run():
    print("\n" + "=" * 70)
    print("AR(1) Demo  |  Part 1: Single-run")
    print(f"n={PARAMS['n']}  p={PARAMS['p']}  s={PARAMS['s']}  rho={PARAMS['rho']:.2f}"
          f"  tFDR={PARAMS['tFDR']:.2f}  amplitude={PARAMS['amplitude']:.2f}")
    print("=" * 70 + "\n")

    support = make_support_random(PARAMS["s"], PARAMS["p"], PARAMS["seed"])
    dat = dgp_ar1_snr(PARAMS["n"], PARAMS["p"], support, amplitude=PARAMS["amplitude"],
                      rho=PARAMS["rho"], snr=PARAMS["snr"], seed=PARAMS["seed"])
    sel = build_da_selector("AR1", dat["X"], dat["y"], PARAMS["tFDR"], PARAMS["K"],
                            "TLARS", **_DA_EXTRA)

    selected = sorted(list(sel.selected_indices))
    tsup = dat["true_support"].tolist()
    n_tp = len(set(selected) & set(tsup))

    print("-" * 70)
    print(f"  True support (0-based): {{{', '.join(map(str, tsup))}}}")
    print(f"  T_stop = {sel.T_stop},  L = {sel.L}")
    print(f"  Selected: {len(selected)}  |  TP = {n_tp}  FP = {len(selected) - n_tp}")
    print(f"  TPP = {compute_tpp(selected, tsup):.3f}  |  "
          f"FDP = {compute_fdp(selected, tsup):.3f}  (target tFDR <= {PARAMS['tFDR']:.2f})")
    print("-" * 70 + "\n")


# ==============================================================================
# Part 2: MC SNR sweep (CappedSpread + Random support)
# ==============================================================================

def part_2_snr_sweep():
    snr_grid = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

    def sweep(support, label):
        print("=" * 70)
        print(f"AR(1) MC — SNR sweep  |  support: {label}")
        print(f"n={MC['n']}  p={MC['p']}  s={MC['s']}  rho={MC['rho']:.1f}  {MC['num_MC']} MC")
        print("=" * 70 + "\n")
        rows = []
        for snr in snr_grid:
            print(f"  [SNR = {snr:.2f}] running {MC['num_MC']} MC trials ...")
            r = run_mc(MC["n"], MC["p"], support, MC["amplitude"], snr, MC["tFDR"],
                       MC["K"], MC["num_MC"], MC["seed"], "dgp_ar1_snr",
                       {"rho": MC["rho"]}, "AR1", _DA_EXTRA, s=MC["s"])
            rows.append(dict(SNR=snr, **r))
        return rows

    capped = make_support_capped_spread(MC["s"], MC["p"], MC["max_gap"])
    print(f"\n  Results — CappedSpread(max_gap={MC['max_gap']}):")
    print_table(sweep(capped, f"CappedSpread(max_gap={MC['max_gap']})"), "SNR")
    print("\n  Results — Random support (redrawn per trial):")
    print_table(sweep(None, "Random (per-trial)"), "SNR")


# ==============================================================================
# Part 3: MC rho sweep (CappedSpread + Random support)
# ==============================================================================

def part_3_rho_sweep():
    rho_grid = [round(0.1 * i, 1) for i in range(10)]   # 0.0 .. 0.9

    def sweep(support, label):
        print("=" * 70)
        print(f"AR(1) MC — Rho sweep  |  support: {label}")
        print(f"n={MC['n']}  p={MC['p']}  s={MC['s']}  SNR={MC['snr']:.1f}  {MC['num_MC']} MC")
        print("=" * 70 + "\n")
        rows = []
        for rho in rho_grid:
            print(f"  [rho = {rho:.1f}] running {MC['num_MC']} MC trials ...")
            r = run_mc(MC["n"], MC["p"], support, MC["amplitude"], MC["snr"], MC["tFDR"],
                       MC["K"], MC["num_MC"], MC["seed"], "dgp_ar1_snr",
                       {"rho": rho}, "AR1", _DA_EXTRA, s=MC["s"])
            rows.append(dict(rho=rho, **r))
        return rows

    capped = make_support_capped_spread(MC["s"], MC["p"], MC["max_gap"])
    print(f"\n  Results — CappedSpread(max_gap={MC['max_gap']}):")
    print_table(sweep(capped, f"CappedSpread(max_gap={MC['max_gap']})"), "rho",
                value_fmt="{:<6.1f}", width=6)
    print("\n  Results — Random support (redrawn per trial):")
    print_table(sweep(None, "Random (per-trial)"), "rho", value_fmt="{:<6.1f}", width=6)


# ==============================================================================
# Part 4: 2D Gap x Rho sweep (kappa boundary)
# ==============================================================================

def part_4_gap_rho_sweep():
    gap_grid = [100, 50, 20, 15, 10, 5, 1]
    rho_grid = [round(0.1 * i, 1) for i in range(10)]
    with np.errstate(divide="ignore"):
        kappa = np.ceil(np.log(0.02) / np.log(np.asarray(rho_grid))).astype(int)
    kappa[np.asarray(rho_grid) == 0.0] = 0

    print("=" * 70)
    print("AR(1) Demo  |  Part 4: Gap x Rho sweep")
    print(f"n={MC['n']}  p={MC['p']}  s={MC['s']}  SNR={MC['snr']:.1f}  {MC['num_MC']} MC")
    print(f"  gap_grid: {{{', '.join(map(str, gap_grid))}}}")
    print("=" * 70 + "\n")

    tpp_cs = np.full((len(rho_grid), len(gap_grid)), np.nan)
    fdp_cs = np.full_like(tpp_cs, np.nan)
    rand_tpp = np.zeros(len(rho_grid))
    rand_fdp = np.zeros(len(rho_grid))

    print("  [A] CappedSpread(gap) sweep ...\n")
    for i, rho in enumerate(rho_grid):
        for j, gap in enumerate(gap_grid):
            support = make_support_capped_spread(MC["s"], MC["p"], gap)
            print(f"  rho={rho:.1f} gap={gap:3d} kappa={kappa[i]:<3d} running {MC['num_MC']} MC ...")
            r = run_mc(MC["n"], MC["p"], support, MC["amplitude"], MC["snr"], MC["tFDR"],
                       MC["K"], MC["num_MC"], MC["seed"], "dgp_ar1_snr",
                       {"rho": rho}, "AR1", _DA_EXTRA, s=MC["s"])
            tpp_cs[i, j] = r["mean_TPP"]
            fdp_cs[i, j] = r["mean_FDP"]

    print("\n  [B] Random support sweep (redrawn each trial) ...\n")
    for i, rho in enumerate(rho_grid):
        print(f"  rho={rho:.1f} kappa={kappa[i]:<3d} Random support, running {MC['num_MC']} MC ...")
        r = run_mc(MC["n"], MC["p"], None, MC["amplitude"], MC["snr"], MC["tFDR"],
                   MC["K"], MC["num_MC"], MC["seed"], "dgp_ar1_snr",
                   {"rho": rho}, "AR1", _DA_EXTRA, s=MC["s"])
        rand_tpp[i] = r["mean_TPP"]
        rand_fdp[i] = r["mean_FDP"]

    print("\n  DA+AR1 correction window half-width: kappa = ceil(log(0.02) / log(rho))")
    print("  rho   : " + "  ".join(f"{r:<5.1f}" for r in rho_grid))
    print("  kappa : " + "  ".join(f"{k:<5d}" for k in kappa))

    row_labels = [f"rho={r:.1f}" for r in rho_grid]
    col_labels = [f"gap={g}" for g in gap_grid] + ["Random"]
    tpp_comb = np.column_stack([tpp_cs, rand_tpp])
    fdp_comb = np.column_stack([fdp_cs, rand_fdp])
    print_matrix(tpp_comb, "mean_TPP  (cols: CappedSpread gap values + Random)",
                 row_labels, col_labels)
    print_matrix(fdp_comb, "mean_FDP  (cols: CappedSpread gap values + Random)",
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
        part_4_gap_rho_sweep()
    print("\ndemo_trex_da_01_ar1.py complete.")
