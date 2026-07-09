# ==============================================================================
# demo_trex_da_02_equi.py
# ==============================================================================
#
# DA-TRex demo and MC simulation for the full-matrix equicorrelated DGP.
#
# Reuses the two-level factor-model DGP (dgp_bt_snr) with rho_within ==
# rho_between, which collapses to pure equicorrelation (compound symmetry).
# The selector uses the EQUI dependency-aware method.
#
#   Part 1: single-run demo.
#   Part 2: MC SNR sweep (rho = 0.25).
#   Part 3: MC rho sweep (SNR = 2.0, rho_within == rho_between).
#
# Mirrors R/trex_selector_methods/trex_da/demo_trex_da_02_equi.R and the equi
# parts of cpp/.../demo_trex_da_02_mc_sim_equi_and_bt. Heavy MC parts are opt-in
# via the RUN_PART_* flags (Part 1 runs by default).
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from dgp_generators import dgp_bt_snr
from da_sim_common import (make_support_random, build_da_selector, run_mc,
                           print_table, compute_fdp, compute_tpp)

RUN_PART_1 = True
RUN_PART_2 = False
RUN_PART_3 = False

MC = dict(n=300, p=1000, s=10, amplitude=3.0, n_blocks=10, tFDR=0.2,
          K=20, num_MC=200, seed=2026)


def _run_equi_mc(rho, snr):
    return run_mc(MC["n"], MC["p"], None, MC["amplitude"], snr, MC["tFDR"],
                  MC["K"], MC["num_MC"], MC["seed"], "dgp_bt_snr",
                  {"n_blocks": MC["n_blocks"], "rho_within": rho, "rho_between": rho},
                  "EQUI", {}, s=MC["s"])


# ==============================================================================
# Part 1: single-run demo
# ==============================================================================

def part_1_single_run():
    print("\n" + "=" * 70)
    print("EQUI Demo  |  Part 1: Single-run")
    print("n=150  p=500  s=5  rho=0.25  SNR=2.0")
    print("=" * 70 + "\n")

    support = make_support_random(5, 500, 2026)
    dat = dgp_bt_snr(150, 500, support, amplitude=3.0, n_blocks=5,
                     rho_within=0.25, rho_between=0.25, snr=2.0, seed=2026)
    sel = build_da_selector("EQUI", dat["X"], dat["y"], 0.2, 20, "TLARS")

    selected = sorted(list(sel.selected_indices))
    tsup = dat["true_support"].tolist()
    n_tp = len(set(selected) & set(tsup))

    print("-" * 70)
    print(f"  True support (0-based): {{{', '.join(map(str, tsup))}}}")
    print(f"  Selected: {len(selected)}  |  TP = {n_tp}  FP = {len(selected) - n_tp}")
    print(f"  TPP = {compute_tpp(selected, tsup):.3f}  |  "
          f"FDP = {compute_fdp(selected, tsup):.3f}  (target tFDR <= 0.20)")
    print("-" * 70 + "\n")


# ==============================================================================
# Part 2: MC SNR sweep (rho = 0.25)
# ==============================================================================

def part_2_snr_sweep():
    rho = 0.25
    snr_grid = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
    print("=" * 70)
    print(f"EQUI MC — SNR sweep  |  rho={rho:.2f}  {MC['num_MC']} MC")
    print("=" * 70 + "\n")
    rows = []
    for snr in snr_grid:
        print(f"  [SNR = {snr:.2f}] running {MC['num_MC']} MC trials ...")
        rows.append(dict(SNR=snr, **_run_equi_mc(rho, snr)))
    print("\n  Results — SNR sweep:")
    print_table(rows, "SNR", value_fmt="{:<8.2f}", width=8)


# ==============================================================================
# Part 3: MC rho sweep (SNR = 2.0, rho_within == rho_between)
# ==============================================================================

def part_3_rho_sweep():
    snr = 2.0
    rho_grid = [round(0.1 * i, 1) for i in range(10)]
    print("=" * 70)
    print(f"EQUI MC — 1D Pure Equi Sweep  |  SNR={snr:.1f}  {MC['num_MC']} MC")
    print("  Sweeping rho (rho_within == rho_between)")
    print("=" * 70 + "\n")
    rows = []
    for rho in rho_grid:
        print(f"  [rho = {rho:.1f}] running {MC['num_MC']} MC trials ...")
        rows.append(dict(rho=rho, **_run_equi_mc(rho, snr)))
    print("\n  Results — 1D Pure Equi Sweep:")
    print_table(rows, "rho", value_fmt="{:<8.2f}", width=8)


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
    print("\ndemo_trex_da_02_equi.py complete.")
