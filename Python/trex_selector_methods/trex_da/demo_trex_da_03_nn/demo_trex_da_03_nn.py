# ==============================================================================
# demo_trex_da_03_nn.py
# ==============================================================================
#
# DA-TRex demo and MC simulation for the Nearest-Neighbours (NN) / MA(kappa)
# banded DGP: Sigma[j, k] != 0 only for |j - k| <= kappa.
#
#   Part 1: single-run demo.
#   Part 2: SNR sweep.
#   Part 3: rho sweep.
#   Part 4: kappa sweep.
#   Part 5: 2D kappa x rho sweep.
#   Part 6: 2D SNR x rho sweep.
#
# Mirrors R/trex_selector_methods/trex_da/demo_trex_da_04_nn_01_ma.R and
# cpp/.../demo_trex_da_03_mc_sim_nn. Parts 1 and 6 are active by default (as in
# R); the other sweeps are opt-in via the RUN_PART_* flags.
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

from dgp_generators import dgp_nn_snr
from da_sim_common import (make_support_capped_spread, make_support_random,
                           build_da_selector, run_mc, print_table, print_matrix,
                           compute_fdp, compute_tpp)

RUN_PART_1 = True
RUN_PART_2 = False
RUN_PART_3 = False
RUN_PART_4 = False
RUN_PART_5 = False
RUN_PART_6 = True

MC = dict(n=300, p=1000, s=10, amplitude=3.0, kappa=3, rho=0.7, snr=2.0,
          tFDR=0.2, K=20, num_MC=200, seed=2026, max_gap=20)


def _run_nn_mc(support, kappa, rho, snr):
    return run_mc(MC["n"], MC["p"], support, MC["amplitude"], snr, MC["tFDR"],
                  MC["K"], MC["num_MC"], MC["seed"], "dgp_nn_snr",
                  {"rho": rho, "kappa": kappa}, "NN", {}, s=MC["s"])


# ==============================================================================
# Part 1: single-run demo
# ==============================================================================

def part_1_single_run():
    n, p, s, kappa, rho, snr, tFDR, K = 150, 500, 5, 1, 0.7, 2.0, 0.2, 20
    print("\n" + "=" * 70)
    print("NN Demo  |  Part 1: Single-run")
    print(f"n={n}  p={p}  s={s}  kappa={kappa}  rho={rho:.1f}  SNR={snr:.1f}")
    print("=" * 70 + "\n")

    support = make_support_random(s, p, 2026)
    dat = dgp_nn_snr(n, p, support, amplitude=3.0, kappa=kappa, rho=rho, snr=snr, seed=2026)
    sel = build_da_selector("NN", dat["X"], dat["y"], tFDR, K, "TLARS")

    selected = sorted(list(sel.selected_indices))
    tsup = dat["true_support"].tolist()
    n_tp = len(set(selected) & set(tsup))
    print("-" * 70)
    print(f"  True support (0-based): {{{', '.join(map(str, tsup))}}}")
    print(f"  T_stop = {sel.T_stop},  L = {sel.L}")
    print(f"  Selected: {len(selected)}  |  TP = {n_tp}  FP = {len(selected) - n_tp}")
    print(f"  TPP = {compute_tpp(selected, tsup):.3f}  |  "
          f"FDP = {compute_fdp(selected, tsup):.3f}  (target tFDR <= {tFDR:.2f})")
    print("-" * 70 + "\n")


# ==============================================================================
# Parts 2-4: 1D sweeps
# ==============================================================================

def _sweep_1d(title, param_name, grid, value_fmt, width, make_kwargs):
    print("=" * 70)
    print(title + f"  |  {MC['num_MC']} MC")
    print("=" * 70 + "\n")
    support = make_support_capped_spread(MC["s"], MC["p"], MC["max_gap"])
    rows = []
    for v in grid:
        print(f"  [{param_name} = {v}] running {MC['num_MC']} MC trials ...")
        rows.append({param_name: v, **_run_nn_mc(support, **make_kwargs(v))})
    print_table(rows, param_name, value_fmt=value_fmt, width=width)


def part_2_snr_sweep():
    _sweep_1d("NN MC — SNR sweep", "SNR", [0.1, 0.2, 0.5, 1.0, 2.0, 5.0],
              "{:<8.2f}", 8, lambda snr: dict(kappa=MC["kappa"], rho=MC["rho"], snr=snr))


def part_3_rho_sweep():
    _sweep_1d("NN MC — rho sweep", "rho", [round(0.1 * i, 1) for i in range(1, 10)],
              "{:<6.1f}", 6, lambda rho: dict(kappa=MC["kappa"], rho=rho, snr=MC["snr"]))


def part_4_kappa_sweep():
    _sweep_1d("NN MC — kappa sweep", "kappa", [1, 2, 3, 5, 7, 10, 15],
              "{:<6d}", 6, lambda k: dict(kappa=k, rho=MC["rho"], snr=MC["snr"]))


# ==============================================================================
# Parts 5-6: 2D sweeps
# ==============================================================================

def _sweep_2d(title, row_name, row_grid, col_name, col_grid, cell_kwargs):
    print("=" * 70)
    print(title + f"  |  {MC['num_MC']} MC")
    print("=" * 70 + "\n")
    support = make_support_capped_spread(MC["s"], MC["p"], MC["max_gap"])
    tpp = np.zeros((len(row_grid), len(col_grid)))
    fdp = np.zeros_like(tpp)
    for i, rv in enumerate(row_grid):
        for j, cv in enumerate(col_grid):
            print(f"  [{row_name}={rv}  {col_name}={cv}] running {MC['num_MC']} MC ...")
            r = _run_nn_mc(support, **cell_kwargs(rv, cv))
            tpp[i, j] = r["mean_TPP"]
            fdp[i, j] = r["mean_FDP"]
    rl = [f"{row_name}={rv}" for rv in row_grid]
    cl = [f"{col_name}={cv}" for cv in col_grid]
    print_matrix(tpp, "mean_TPP", rl, cl)
    print_matrix(fdp, "mean_FDP", rl, cl)


def part_5_kappa_rho_sweep():
    _sweep_2d("NN MC — 2D kappa x rho sweep (SNR=2.0)",
              "kappa", [1, 2, 3, 5, 7, 10, 15], "rho", [0.1, 0.3, 0.5, 0.7, 0.8, 0.9],
              lambda k, rho: dict(kappa=k, rho=rho, snr=2.0))


def part_6_snr_rho_sweep():
    _sweep_2d("NN MC — 2D SNR x rho sweep (kappa=3)",
              "SNR", [0.1, 0.5, 1.0, 2.0, 5.0],
              "rho", [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9],
              lambda snr, rho: dict(kappa=MC["kappa"], rho=rho, snr=snr))


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
        part_4_kappa_sweep()
    if RUN_PART_5:
        part_5_kappa_rho_sweep()
    if RUN_PART_6:
        part_6_snr_rho_sweep()
    print("\ndemo_trex_da_03_nn.py complete.")
