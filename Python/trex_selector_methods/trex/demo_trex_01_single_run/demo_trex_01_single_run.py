# ==============================================================================
# demo_trex_01_single_run.py
# ==============================================================================
#
# Classical T-Rex Selector demo — single-run.
# Mirrors R/trex_selector_methods/trex/demo_trex_01_single_run.R and
# cpp/trex_selector_methods/trex/demo_trex_01_single_run.cpp.
#
# Demo content:
#  - Low- and high-dimensional settings.
#  - Fixed true support; TLARS solver.
#  - Reports phi_prime, phi_mat, fdp_hat_mat, r_mat, voting_grid.
#
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import numpy as np
import trex_selector_neo as ts

from dgp_gauss_snr import dgp_gauss_snr
from trex_sim_common import _make_trex_ctrl_from_dict
from trex_helpers import print_single_run_result

# ==============================================================================
# Global parameters
# ==============================================================================

# Fixed true support (0-based, same as C++ demo).
# C++ (0-based): {27, 149, 43, 128, 42, 4}
# R   (1-based): {28, 150, 44, 129, 43, 5}
TRUE_SUPPORT = np.array([27, 149, 43, 128, 42, 4], dtype=np.intp)


# ==============================================================================
# Part 1: single-run basic demo
# ==============================================================================

def part_1_single_run():
    tFDR = 0.1
    snr  = 1.0
    seed = 4231

    ctrl_dict = dict(
        solver_name               = "TLARS",
        K                         = 20,
        max_dummy_multiplier      = 10,
        use_max_T_stop            = True,
        lloop_strategy            = "HCONCAT",
        tloop_stagnation_stop     = False,
        tloop_max_stagnant_steps  = 5,
        parallel_rnd_experiments  = False,   # mirrors R use_openmp=FALSE
        lambda2                   = 0.1,
        rho_afs                   = 0.3,
        ncgmp_variant             = 0,
    )

    def run_single(high_dim):
        n         = 150  if high_dim else 5000
        p         = 300  if high_dim else 1000
        dim_label = "High-dimensional (p > n)" if high_dim else "Low-dimensional (n > p)"

        print(f"\n{'=' * 70}")
        print(f"Part 1: Single-run basic T-Rex demo  |  {dim_label}")
        print(f"n = {n},  p = {p},  tFDR = {tFDR:.2f},  SNR = {snr:.1f}\n")

        print("Generating synthetic data ...")
        dat = dgp_gauss_snr(n, p, TRUE_SUPPORT, amplitude=1.0, snr=snr, seed=seed)

        print("Running T-Rex Selector (TLARS) ...\n")
        ctrl = _make_trex_ctrl_from_dict(ctrl_dict)
        sel  = ts.TRexSelector(
            np.asfortranarray(dat["X"]),
            dat["y"],
            tFDR=tFDR,
            trex_control=ctrl,
            seed=seed,
            verbose=True,
        )
        sel.select()

        print_single_run_result(
            f"Part 1 \u2014 i.i.d. Gaussian  [{dim_label}, TLARS]",
            dat, sel, tFDR,
        )

    run_single(high_dim=False)
    run_single(high_dim=True)


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    part_1_single_run()
    print("\ndemo_trex_01_single_run.py complete.")
