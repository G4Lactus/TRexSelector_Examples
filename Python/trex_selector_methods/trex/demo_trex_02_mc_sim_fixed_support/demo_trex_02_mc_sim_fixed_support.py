# ==============================================================================
# demo_trex_02_mc_sim_fixed_support.py
# ==============================================================================
#
# Classical T-Rex selector demo — Monte Carlo simulation with fixed true support.
# Mirrors R/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support.R and
# cpp/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support.cpp.
#
# Demo content:
#  - 14 solvers; SNR sweep {0.1, 0.5, 1.0, 2.0, 5.0}.
#  - Fixed support drawn once with seed 24 (shared across all trials).
#  - Reports averaged FDR / TPR per solver x SNR.
#
# Python downscale: num_MC = 100 (C++ runs 200) to keep the
# ProcessPoolExecutor-based run time moderate.
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

from support_generators import make_support_random
from trex_sim_common import (
    SOLVERS_DEFAULT,
    run_mc_trex,
    print_solver_table,
    save_mc_results,
)

# ==============================================================================
# Global parameters
# ==============================================================================

_OUT_DIR     = os.path.join(_THIS_DIR, "simulation_results", "data")
_NUM_WORKERS = 6
_NUM_MC      = 100


# ==============================================================================
# MC simulation — fixed true support, solver comparison
# ==============================================================================

def trx_02_mc_sim_fixed_support(num_MC=_NUM_MC, num_workers=_NUM_WORKERS, high_dim=True):

    n = 300  if high_dim else 1000
    p = 1000 if high_dim else 300
    s = 10
    tFDR        = 0.1
    snr_values  = [0.1, 0.5, 1.0, 2.0, 5.0]
    stagnant_sw = 7

    dim_label = "High-dimensional (p > n)" if high_dim else "Low-dimensional (n > p)"
    print(f"\n{'=' * 70}")
    print(f"MC Simulation \u2014 Fixed True Support  |  {dim_label}")
    print(f"n = {n},  p = {p},  s = {s},  num_MC = {num_MC}")

    # Base control shared by all solvers
    base_ctrl = dict(
        K                        = 20,
        max_dummy_multiplier     = 10,
        use_max_T_stop           = True,
        lloop_strategy           = "STANDARD",
        tloop_stagnation_stop    = True,
        tloop_max_stagnant_steps = stagnant_sw,
        opt_threshold            = 0.75,
        parallel_rnd_experiments = False,   # mirrors R use_openmp=FALSE
        use_memory_mapping       = False,
    )

    # Fixed support drawn once — shared across all trials and solvers
    true_support = make_support_random(s, p, seed=24)
    print(f"True support (cardinality {s}): {{{', '.join(str(x) for x in true_support)}}}\n")

    n_sol = len(SOLVERS_DEFAULT)
    n_snr = len(snr_values)
    fdr_mat = np.zeros((n_sol, n_snr))
    tpr_mat = np.zeros((n_sol, n_snr))

    # Simulation loop: iterate over solvers, then SNR values
    for si, solver in enumerate(SOLVERS_DEFAULT):
        print(f"Solver: {solver['name']}")
        for snr_idx, snr in enumerate(snr_values):
            base_seed = 24 + snr_idx * 1000   # mirrors R: 24L + (snr_idx-1L)*1000L (1-based)
            lbl = f"SNR={snr:.2f} [{solver['name']}]"
            result = run_mc_trex(
                n=n, p=p, s=s, snr=snr,
                solver_desc=solver,
                base_ctrl_dict=base_ctrl,
                tFDR=tFDR,
                base_seed=base_seed,
                num_MC=num_MC,
                max_workers=num_workers,
                fixed_support=true_support,
                track_L_T=False,
                label=lbl,
            )
            fdr_mat[si, snr_idx] = result["mean_FDR"]
            tpr_mat[si, snr_idx] = result["mean_TPR"]
        print()

    solver_names = [s["name"] for s in SOLVERS_DEFAULT]
    stem = (f"demo_trex_02_mc_sim_fixed_support_trex_results_n{n}_p{p}"
            f"_stagnation_window_{stagnant_sw}")
    save_mc_results(
        _OUT_DIR, stem, num_MC,
        solver_names, snr_values,
        fdr_mat, tpr_mat,
    )


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    trx_02_mc_sim_fixed_support()
    print("\ndemo_trex_02_mc_sim_fixed_support.py complete.")
