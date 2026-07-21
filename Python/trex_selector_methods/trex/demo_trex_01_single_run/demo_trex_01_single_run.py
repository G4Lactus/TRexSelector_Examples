# ==============================================================================
# demo_trex_01_single_run.py
# ==============================================================================
#
# Classical T-Rex Selector demo — single-run.
# Mirrors cpp/trex_selector_methods/trex/demo_trex_01_single_run.cpp and
# R/trex_selector_methods/trex/demo_trex_01_single_run.R.
#
# Demo content:
#  - High- and low-dimensional settings (run in that order, as in C++).
#  - Fixed true support; fixed +1 coefficients (rnd_coef=False); TLARS solver.
#  - Reports selected indices, FDP/TPP, Phi_prime and Phi (printed sparsely),
#    FDP_hat, R, and the voting grid — dual console + file output to
#    simulation_results/data/d01_trex_basic_n{n}_p{p}.txt.
#
# Seed note: data seed 58 and selector seed 42 mirror the C++ demo
# structurally; the RNG streams differ across languages, so the drawn data
# (and hence the exact selection) are not bit-identical to the C++ output.
#
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

import trex_selector_neo as tsn

from dgp_gauss_snr import dgp_gauss_snr
from trex_helpers import make_dual_out, print_single_run_result
from trex_sim_common import _make_trex_ctrl_from_dict

# ==============================================================================
# Global parameters
# ==============================================================================

# Fixed true support (0-based, same as C++ demo).
# C++ (0-based): {27, 149, 43, 128, 42, 4}
# R   (1-based): {28, 150, 44, 129, 43, 5}
TRUE_SUPPORT = np.array([27, 149, 43, 128, 42, 4], dtype=np.intp)

_OUT_DIR = os.path.join(_THIS_DIR, "simulation_results", "data")


# ==============================================================================
# Part 1: single-run basic demo
# ==============================================================================

def part_1_single_run():
    tFDR      = 0.1
    snr       = 1.0
    data_seed = 58   # mirrors C++ SyntheticData(..., seed=58)
    trex_seed = 42   # mirrors C++ TRexSelector(..., 42, true)

    ctrl_dict = dict(
        solver_name               = "TLARS",
        K                         = 20,
        max_dummy_multiplier      = 10,
        use_max_T_stop            = True,
        lloop_strategy            = "STANDARD",
        tloop_stagnation_stop     = False,
        tloop_max_stagnant_steps  = 5,
        parallel_rnd_experiments  = False,   # mirrors C++ / R use_openmp=FALSE
        lambda2                   = 0.1,
        rho_afs                   = 0.3,
        ncgmp_variant             = 0,
    )

    def run_single(high_dim, rnd_coef=False):
        n         = 150 if high_dim else 5000
        p         = 300 if high_dim else 1000
        dim_label = "High-dimensional (p > n)" if high_dim else "Low-dimensional (n > p)"

        # The C++ rnd_coef path selects a fixed hardcoded vector that mimics
        # randomly drawn coefficients (it is NOT a random draw); main() runs
        # with rnd_coef=False, i.e. all coefficients equal 1.
        coefs = (np.array([-0.4, -0.25, -0.8, 1.1, 2.5, -1.2])
                 if rnd_coef else np.ones(len(TRUE_SUPPORT)))

        os.makedirs(_OUT_DIR, exist_ok=True)
        out_path = os.path.join(_OUT_DIR, f"d01_trex_basic_n{n}_p{p}.txt")

        with open(out_path, "w") as fh:
            out = make_dual_out(fh)

            out(f"\n{'=' * 70}\n")
            out(f"Part 1: Single-run basic T-Rex demo  |  {dim_label}\n")
            out(f"n = {n},  p = {p},  tFDR = {tFDR:.2f},  SNR = {snr:.1f}\n\n")

            out("Generating synthetic data ...\n")
            dat = dgp_gauss_snr(n, p, TRUE_SUPPORT, coefs=coefs, snr=snr,
                                seed=data_seed)

            out("Running T-Rex Selector (TLARS) ...\n\n")
            ctrl = _make_trex_ctrl_from_dict(ctrl_dict)
            sel  = tsn.TRexSelector(
                np.asfortranarray(dat["X"]),
                dat["y"],
                tFDR=tFDR,
                trex_control=ctrl,
                seed=trex_seed,
                verbose=True,
            )
            sel.select()

            print_single_run_result(
                f"Part 1 — i.i.d. Gaussian  [{dim_label}, TLARS]",
                dat, sel, tFDR, out=out,
            )

        print(f"[Info] Results saved to: {out_path}")

    # C++ main() runs the high-dimensional scenario first, then the
    # low-dimensional one; both with fixed +1 coefficients.
    run_single(high_dim=True)
    run_single(high_dim=False)


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    part_1_single_run()
    print("\ndemo_trex_01_single_run.py complete.")
