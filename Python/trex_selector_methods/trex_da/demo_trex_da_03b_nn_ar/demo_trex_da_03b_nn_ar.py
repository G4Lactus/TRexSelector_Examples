# ==============================================================================
# demo_trex_da_03b_nn_ar.py
# ==============================================================================
#
# TREX+DA+NN applied to AR(1) data (method-mismatch test): the NN / MA-band
# correction is run on AR(1) Toeplitz designs to probe how it behaves when the
# true dependence is AR(1) rather than MA(kappa).
#
#   Part 1: single-run demo.
#   Part 2: SNR sweep.
#   Part 3: rho sweep.
#   Part 4: 2D SNR x rho sweep.
#
# Mirrors R/trex_selector_methods/trex_da/demo_trex_da_04_nn_02_ar.R and
# cpp/.../demo_trex_da_03b_mc_sim_nn_ar. Part 4 is active by default (as in R);
# the rest are opt-in via the RUN_PART_* flags.
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

RUN_PART_1 = False
RUN_PART_2 = False
RUN_PART_3 = False
RUN_PART_4 = True

MC = dict(n=300, p=1000, s=10, amplitude=3.0, rho=0.7, snr=2.0, tFDR=0.2,
          K=20, num_MC=200, seed=2026, max_gap=20)


def _run_nn_on_ar(support, rho, snr):
    return run_mc(MC["n"], MC["p"], support, MC["amplitude"], snr, MC["tFDR"],
                  MC["K"], MC["num_MC"], MC["seed"], "dgp_ar1_snr",
                  {"rho": rho}, "NN", {}, s=MC["s"])


def part_1_single_run():
    n, p, s, rho, snr, tFDR, K = 150, 500, 5, 0.7, 2.0, 0.2, 20
    print("\n" + "=" * 70)
    print("AR(1) Demo | Part 1: Single-run  [trex+DA+NN on AR(1) data]")
    print(f"n={n}  p={p}  s={s}  rho={rho:.1f}  SNR={snr:.1f}")
    print("=" * 70 + "\n")

    support = make_support_random(s, p, 2026)
    dat = dgp_ar1_snr(n, p, support, amplitude=3.0, rho=rho, snr=snr, seed=2026)
    sel = build_da_selector("NN", dat["X"], dat["y"], tFDR, K, "TLARS")

    selected = sorted(list(sel.selected_indices))
    tsup = dat["true_support"].tolist()
    n_tp = len(set(selected) & set(tsup))
    print("-" * 70)
    print(f"  T_stop = {sel.T_stop},  L = {sel.L}")
    print(f"  Selected: {len(selected)}  |  TP = {n_tp}  FP = {len(selected) - n_tp}")
    print(f"  TPP = {compute_tpp(selected, tsup):.3f}  |  "
          f"FDP = {compute_fdp(selected, tsup):.3f}  (target tFDR <= {tFDR:.2f})")
    print("-" * 70 + "\n")


def _snr_sweep(support, label):
    print("=" * 70)
    print(f"AR(1) MC — SNR sweep (NN method)  |  support: {label}  {MC['num_MC']} MC")
    print("=" * 70 + "\n")
    rows = []
    for snr in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        print(f"  [SNR = {snr:.2f}] running {MC['num_MC']} MC trials ...")
        rows.append(dict(SNR=snr, **_run_nn_on_ar(support, MC["rho"], snr)))
    print_table(rows, "SNR", value_fmt="{:<8.2f}", width=8)


def part_2_snr_sweep():
    _snr_sweep(make_support_capped_spread(MC["s"], MC["p"], MC["max_gap"]),
               f"CappedSpread(max_gap={MC['max_gap']})")


def part_3_rho_sweep():
    support = make_support_capped_spread(MC["s"], MC["p"], MC["max_gap"])
    print("=" * 70)
    print(f"AR(1) MC — rho sweep (NN method)  |  {MC['num_MC']} MC")
    print("=" * 70 + "\n")
    rows = []
    for rho in [round(0.1 * i, 1) for i in range(1, 10)]:
        print(f"  [rho = {rho:.1f}] running {MC['num_MC']} MC trials ...")
        rows.append(dict(rho=rho, **_run_nn_on_ar(support, rho, MC["snr"])))
    print_table(rows, "rho", value_fmt="{:<6.1f}", width=6)


def part_4_snr_rho_sweep():
    snr_grid = [2.0, 5.0]
    rho_grid = [round(0.1 * i, 1) for i in range(1, 10)]
    print("=" * 70)
    print(f"AR(1) MC — 2D SNR x Rho sweep (NN method)  |  {MC['num_MC']} MC")
    print("=" * 70 + "\n")
    support = make_support_capped_spread(MC["s"], MC["p"], MC["max_gap"])
    tpp = np.zeros((len(snr_grid), len(rho_grid)))
    fdp = np.zeros_like(tpp)
    for i, snr in enumerate(snr_grid):
        for j, rho in enumerate(rho_grid):
            print(f"  [SNR={snr}  rho={rho}] running {MC['num_MC']} MC ...")
            r = _run_nn_on_ar(support, rho, snr)
            tpp[i, j] = r["mean_TPP"]
            fdp[i, j] = r["mean_FDP"]
    rl = [f"SNR={snr}" for snr in snr_grid]
    cl = [f"rho={rho}" for rho in rho_grid]
    print_matrix(tpp, "mean_TPP", rl, cl)
    print_matrix(fdp, "mean_FDP", rl, cl)


if __name__ == "__main__":
    if RUN_PART_1:
        part_1_single_run()
    if RUN_PART_2:
        part_2_snr_sweep()
    if RUN_PART_3:
        part_3_rho_sweep()
    if RUN_PART_4:
        part_4_snr_rho_sweep()
    print("\ndemo_trex_da_03b_nn_ar.py complete.")
