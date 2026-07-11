# ==============================================================================
# demo_trex_scr_03_screen_correlated.py
# ==============================================================================
#
# Screen-TRex Selector, part 03: correlated designs & dependency-aware screening.
#
# Ordinary screening thresholds each variable's voting proportion Phi_j > 0.5
# independently. Under correlated predictors, redundant correlated variables
# *split* the voting evidence, so Phi drops below 0.5 for all of them and power
# collapses. The dependency-aware (DA) variants pre-group strongly correlated
# variables before voting so the evidence is not split:
#   * TREX_DA_AR1        - for AR(1) (banded) correlation,
#   * TREX_DA_EQUI       - for equicorrelation (single shared factor),
#   * TREX_DA_BLOCK_EQUI - for block-equicorrelation (per-block factors).
# Each DA variant also reports an Estimated rho (estimated_correlation).
#
# Mirrors:
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_03_screen_trex_correlated/
#       demo_trex_scr_03_screen_trex_correlated.cpp
#   R/trex_selector_methods/trex_screening/demo_trex_scr_03_screen_correlated/
#       demo_trex_scr_03_screen_correlated.R
#
# Demo content:
#   Part 1: Ordinary Screen-TRex on AR(1) data (baseline, no DA correction).
#   Part 2: DA-AR1 Screen-TRex on AR(1) data.
#   Part 3: DA-EQUI Screen-TRex on equicorrelated data.
#   Part 4: DA-BLOCK-EQUI Screen-TRex on block-equicorrelated data.
#   Part 5: MC SNR sweep on AR(1): Ordinary / Bootstrap / +DA-AR1 variants.
#   Part 6: MC SNR sweep on equicorrelated: Ordinary / Bootstrap / +DA-EQUI.
#   Part 7: MC rho sweep on AR(1) (fixed SNR): the four methods above.
#   Part 8: MC rho sweep on equicorrelated (fixed SNR).
#   Part 9: MC rho sweep on block-equicorrelated (fixed SNR).
#
# The C++ demo uses num_MC = 500; here we use a smaller count. Support indices
# are 0-based.
# ==============================================================================

import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import trex_selector_neo as ts
from trex_scr_common import (dgp_ar1, dgp_equi, dgp_block_equi, make_scr_trex_ctrl,
                             make_screen_ctrl, print_scr_result, run_mc_screen,
                             save_and_print_scr_mc)

OUT_DIR = os.path.join(_THIS_DIR, "simulation_results", "data")

RUN_SINGLE_PARTS = True   # Parts 1-4
RUN_MC_PARTS = True       # Parts 5-9

N, P, P1 = 300, 1000, 10
SNR_GRID = [0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0]
RHO_GRID = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
NUM_MC = 10   # C++ uses 500; reduced so all 5 MC parts finish quickly.


def run_single(title, dat, trex_method, bootstrap=False):
    """Single-run helper: run Screen-TRex, print block."""
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70 + "\n")
    sc = make_screen_ctrl(trex_method=trex_method, bootstrap=bootstrap)
    sel = ts.TRexScreeningSelector(dat["X"], dat["y"], screen_control=sc,
                                   trex_control=make_scr_trex_ctrl(),
                                   seed=42, verbose=False)
    res = sel.select()
    print_scr_result(sel, res, dat["true_support"])


def run_mc_sweep(part_title, dgp_fn, sweep_name, sweep_values, methods, file_stem):
    """Generic MC sweep over one variable (SNR or rho), for a list of methods."""
    print("\n" + "=" * 70)
    print(part_title)
    print("=" * 70 + "\n")
    method_names = [m["name"] for m in methods]

    def init():
        return {nm: [0.0] * len(sweep_values) for nm in method_names}
    fdr_map, tpr_map, est_fdr_map = init(), init(), init()
    trex_ctrl = make_scr_trex_ctrl()

    for m in methods:
        print("=" * 70)
        print(f"Method: {m['name']}")
        print("=" * 70)
        for vi, v in enumerate(sweep_values):
            def make_data(seed, v=v):
                return dgp_fn(v, seed)   # v = SNR or rho
            res = run_mc_screen(
                NUM_MC, f"{m['name']}  {sweep_name}={v:.2f}", make_data, trex_ctrl,
                screen_args={"trex_method": m["trex_method"],
                             "bootstrap": m["bootstrap"]},
                base_seed=42 + vi * 1000)
            fdr_map[m["name"]][vi] = res["avg_fdr"]
            tpr_map[m["name"]][vi] = res["avg_tpr"]
            est_fdr_map[m["name"]][vi] = res["avg_est_fdr"]

    save_and_print_scr_mc(NUM_MC, OUT_DIR, file_stem, sweep_values,
                          method_names, fdr_map, tpr_map, est_fdr_map,
                          sweep_label=sweep_name)


# ==============================================================================
# Parts 1-4: single runs on correlated designs
# ==============================================================================

def single_parts():
    run_single("Part 1: Ordinary Screen-TRex on AR(1) data (no DA correction)",
               dgp_ar1(N, P, P1, rho=0.7, snr=5.0, beta_val=3.0, seed=1),
               "TREX")
    run_single("Part 2: DA-AR1 Screen-TRex on AR(1) data",
               dgp_ar1(N, P, P1, rho=0.5, snr=1.0, beta_val=1.0, seed=1),
               "TREX_DA_AR1")
    run_single("Part 3: DA-EQUI Screen-TRex on equicorrelated data",
               dgp_equi(N, P, P1, rho=0.4, snr=1.0, beta_val=1.0, seed=1),
               "TREX_DA_EQUI")
    run_single("Part 4: DA-BLOCK-EQUI Screen-TRex on block-equicorrelated data",
               dgp_block_equi(N, P, P1, rho=0.5, n_blocks=5, snr=1.0,
                              beta_val=1.0, seed=1),
               "TREX_DA_BLOCK_EQUI")


# ==============================================================================
# Parts 5-9: Monte Carlo sweeps
# ==============================================================================

def mc_parts():
    methods_ar1 = [
        {"name": "Screen-TRex Ordinary",    "trex_method": "TREX",        "bootstrap": False},
        {"name": "Screen-TRex Bootstrap",   "trex_method": "TREX",        "bootstrap": True},
        {"name": "Screen-TRex Ord+DA-AR1",  "trex_method": "TREX_DA_AR1", "bootstrap": False},
        {"name": "Screen-TRex Boot+DA-AR1", "trex_method": "TREX_DA_AR1", "bootstrap": True}]
    methods_equi = [
        {"name": "Screen-TRex Ordinary",     "trex_method": "TREX",         "bootstrap": False},
        {"name": "Screen-TRex Bootstrap",    "trex_method": "TREX",         "bootstrap": True},
        {"name": "Screen-TRex Ord+DA-EQUI",  "trex_method": "TREX_DA_EQUI", "bootstrap": False},
        {"name": "Screen-TRex Boot+DA-EQUI", "trex_method": "TREX_DA_EQUI", "bootstrap": True}]
    methods_block = [
        {"name": "Screen-TRex Ordinary",     "trex_method": "TREX",               "bootstrap": False},
        {"name": "Screen-TRex Bootstrap",    "trex_method": "TREX",               "bootstrap": True},
        {"name": "Screen-TRex Ord+DA-BLOCK", "trex_method": "TREX_DA_BLOCK_EQUI", "bootstrap": False},
        {"name": "Screen-TRex Boot+DA-BLOCK","trex_method": "TREX_DA_BLOCK_EQUI", "bootstrap": True}]

    # Part 5: SNR sweep on AR(1) (rho=0.5, beta=1)
    run_mc_sweep("Part 5: MC SNR sweep on AR(1) data (rho=0.5)",
                 lambda snr, seed: dgp_ar1(N, P, P1, 0.5, snr, 1.0, seed),
                 "SNR", SNR_GRID, methods_ar1, "d03_da_ar1_mc_n300_p1000_rho0.50")

    # Part 6: SNR sweep on equicorrelated (rho=0.4, beta=3)
    run_mc_sweep("Part 6: MC SNR sweep on equicorrelated data (rho=0.4)",
                 lambda snr, seed: dgp_equi(N, P, P1, 0.4, snr, 3.0, seed),
                 "SNR", SNR_GRID, methods_equi, "d03_da_equi_mc_n300_p1000_rho0.40")

    # Part 7: rho sweep on AR(1) (SNR=1, beta=1)
    run_mc_sweep("Part 7: MC rho sweep on AR(1) data (SNR=1)",
                 lambda rho, seed: dgp_ar1(N, P, P1, rho, 1.0, 1.0, seed),
                 "rho", RHO_GRID, methods_ar1, "d03_da_ar1_rho_sweep_n300_p1000_snr1.00")

    # Part 8: rho sweep on equicorrelated (SNR=1, beta=3)
    run_mc_sweep("Part 8: MC rho sweep on equicorrelated data (SNR=1)",
                 lambda rho, seed: dgp_equi(N, P, P1, rho, 1.0, 3.0, seed),
                 "rho", RHO_GRID, methods_equi, "d03_da_equi_rho_sweep_n300_p1000_snr1.00")

    # Part 9: rho sweep on block-equicorrelated (SNR=1, beta=3, 5 blocks)
    run_mc_sweep("Part 9: MC rho sweep on block-equicorrelated data (SNR=1)",
                 lambda rho, seed: dgp_block_equi(N, P, P1, rho, 5, 1.0, 3.0, seed),
                 "rho", RHO_GRID, methods_block,
                 "d03_da_block_equi_rho_sweep_n300_p1000_snr1.00")


# ==============================================================================
# Main
# ==============================================================================

def main():
    if RUN_SINGLE_PARTS:
        single_parts()
    if RUN_MC_PARTS:
        mc_parts()
    print("\ndemo_trex_scr_03_screen_correlated.py complete.")


if __name__ == "__main__":
    main()
