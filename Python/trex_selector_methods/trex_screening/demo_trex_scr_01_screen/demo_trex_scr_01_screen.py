# ==============================================================================
# demo_trex_scr_01_screen.py
# ==============================================================================
#
# Screen-TRex Selector, part 01: basic in-memory variable screening.
#
# Screen-TRex extends the classical T-Rex Selector with a fast variable-
# *screening* step for ultra-high-dimensional data. It thresholds the same
# dummy-based voting proportion Phi_j that T-Rex accumulates internally:
#   * Ordinary    rule: select {j : Phi_j > 0.5}  (a majority vote).
#   * Bootstrap-CI rule: build a bootstrap confidence band around the estimated
#     FDR curve (R_boot replicates) and pick the threshold from that band.
#
# Mirrors:
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_01_screen_trex/
#       demo_trex_scr_01_screen_trex.cpp
#   R/trex_selector_methods/trex_screening/demo_trex_scr_01_screen/
#       demo_trex_scr_01_screen.R
#
# Demo content:
#   Part 1: single-run Ordinary Screen-TRex (high-dimensional, fixed support).
#   Part 2: single-run Bootstrap-CI Screen-TRex (same data).
#   Part 3: Monte Carlo SNR sweep comparing Ordinary vs. Bootstrap-CI screening
#           (FDR, TPR, and the procedure's self-reported Estimated FDR).
#
# The C++ demo uses num_MC = 500; here we use a smaller count (noted below) so
# the demo finishes quickly. Results are statistically comparable, not
# bit-identical (different RNG streams). Support indices are 0-based.
# ==============================================================================

import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import trex_selector_neo as ts
from trex_scr_common import (dgp_iid_snr, make_scr_trex_ctrl, make_screen_ctrl,
                             default_scr_methods, print_scr_result,
                             run_mc_screen, save_and_print_scr_mc)

OUT_DIR = os.path.join(_THIS_DIR, "simulation_results")

RUN_PART_1 = True
RUN_PART_2 = True
RUN_PART_3 = True

# Fixed true support (0-based, matches the C++ demo).
TRUE_SUPPORT = [27, 149, 398, 4, 42]

NUM_MC = 60   # C++ uses 500; reduced for a quick runtime.


# ==============================================================================
# Parts 1 & 2: single-run Ordinary and Bootstrap-CI screening
# ==============================================================================

def run_single(bootstrap):
    n, p, snr = 300, 1000, 5.0
    label = ("Bootstrap-CI (Phi via bootstrap band)" if bootstrap
             else "Ordinary (Phi > 0.5)")

    print("\n" + "=" * 70)
    print(f"Screen-TRex single run  |  {label}")
    print(f"High-dimensional (p > n): n = {n}, p = {p}, "
          f"s = {len(TRUE_SUPPORT)}, SNR = {snr:.1f}\n")

    print("Generating synthetic data ...")
    dat = dgp_iid_snr(n, p, TRUE_SUPPORT, coefs=[1.0] * len(TRUE_SUPPORT),
                      snr=snr, seed=58)

    print("Running Screen-TRex selection ...")
    sc = make_screen_ctrl(trex_method="TREX", bootstrap=bootstrap)
    sel = ts.TRexScreeningSelector(dat["X"], dat["y"], screen_control=sc,
                                   trex_control=make_scr_trex_ctrl(),
                                   seed=42, verbose=False)
    res = sel.select()
    print_scr_result(sel, res, dat["true_support"])


# ==============================================================================
# Part 3: Monte Carlo SNR sweep - Ordinary vs. Bootstrap-CI
# ==============================================================================

def run_mc_sweep():
    print("\n" + "=" * 70)
    print("Part 3: Monte Carlo SNR sweep (Ordinary vs. Bootstrap-CI)")
    print(f"High-dimensional: n = 300, p = 1000, s = 10, num_MC = {NUM_MC}\n")

    import numpy as np
    n, p, support_size = 300, 1000, 10
    snr_values = [0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0]
    methods = default_scr_methods()
    method_names = [m["name"] for m in methods]

    def init():
        return {nm: [0.0] * len(snr_values) for nm in method_names}
    fdr_map, tpr_map, est_fdr_map = init(), init(), init()

    trex_ctrl = make_scr_trex_ctrl()

    for m in methods:
        print("=" * 70)
        print(f"Method: {m['name']}")
        print("=" * 70)
        for si, snr in enumerate(snr_values):
            def make_data(seed, snr=snr):
                rng = np.random.default_rng(seed + 500000)   # isolated support RNG
                sup = rng.choice(p, size=support_size, replace=False)
                return dgp_iid_snr(n, p, sup, [1.0] * support_size, snr, seed)
            base_seed = 24 + si * 1000
            res = run_mc_screen(
                NUM_MC, f"{m['name']}  SNR={snr:.2f}", make_data, trex_ctrl,
                screen_args={"trex_method": m["trex_method"],
                             "bootstrap": m["bootstrap"]},
                base_seed=base_seed)
            fdr_map[m["name"]][si] = res["avg_fdr"]
            tpr_map[m["name"]][si] = res["avg_tpr"]
            est_fdr_map[m["name"]][si] = res["avg_est_fdr"]

    save_and_print_scr_mc(NUM_MC, OUT_DIR, "d01_screen_trex_mc_n300_p1000",
                          snr_values, method_names, fdr_map, tpr_map, est_fdr_map)


# ==============================================================================
# Main
# ==============================================================================

def main():
    if RUN_PART_1:
        run_single(bootstrap=False)
    if RUN_PART_2:
        run_single(bootstrap=True)
    if RUN_PART_3:
        run_mc_sweep()
    print("\ndemo_trex_scr_01_screen.py complete.")


if __name__ == "__main__":
    main()
