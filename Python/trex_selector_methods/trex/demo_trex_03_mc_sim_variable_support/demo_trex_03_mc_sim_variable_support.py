# ==============================================================================
# demo_trex_03_mc_sim_variable_support.py
# ==============================================================================
#
# Classical T-Rex selector demo — Monte Carlo simulation with variable support.
# Mirrors R/trex_selector_methods/trex/demo_trex_03_mc_sim_variable_support.R and
# cpp/trex_selector_methods/trex/demo_trex_03_mc_sim_variable_support.cpp.
#
# Demo content:
#  - 14 solvers; SNR sweep {0.1, 0.2, ..., 2.0, 5.0}  (21 values).
#  - Support and (optionally) coefficients redrawn each trial.
#  - Reports averaged FDR / TPR / Avg L / Avg T per solver x SNR.
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

from trex_sim_common import (
    SOLVERS_DEFAULT,
    run_mc_trex,
    save_mc_results,
)

# ==============================================================================
# Global parameters
# ==============================================================================

_OUT_DIR     = os.path.join(_THIS_DIR, "simulation_results")
_NUM_WORKERS = 6
_NUM_MC      = 100


# ==============================================================================
# MC simulation — variable true support, solver comparison
# ==============================================================================

def trx_03_mc_sim_variable_support(num_MC=_NUM_MC, num_workers=_NUM_WORKERS,
                                   high_dim=True, rnd_coef=False):

    n = 300  if high_dim else 1000
    p = 1000 if high_dim else 300
    s = 10
    tFDR        = 0.1
    snr_values  = list(np.round(np.arange(0.1, 2.1, 0.1), 1)) + [5.0]  # 21 values
    stagnant_sw = 3

    dim_label = "High-dimensional (p > n)" if high_dim else "Low-dimensional (n > p)"
    print(f"\n{'=' * 70}")
    print(f"MC Simulation \u2014 Variable Support  |  {dim_label}")
    print(f"n = {n},  p = {p},  s = {s},  num_MC = {num_MC}\n")

    # Base control shared by all solvers
    base_ctrl = dict(
        K                        = 20,
        max_dummy_multiplier     = 10,
        use_max_T_stop           = True,
        lloop_strategy           = "HCONCAT",
        tloop_stagnation_stop    = True,
        tloop_max_stagnant_steps = stagnant_sw,
        parallel_rnd_experiments = False,   # mirrors R use_openmp=FALSE
    )

    n_sol = len(SOLVERS_DEFAULT)
    n_snr = len(snr_values)
    fdr_mat   = np.zeros((n_sol, n_snr))
    tpr_mat   = np.zeros((n_sol, n_snr))
    avg_L_mat = np.zeros((n_sol, n_snr))
    avg_T_mat = np.zeros((n_sol, n_snr))

    # Simulation loop: iterate over solvers and SNR values
    for si, solver in enumerate(SOLVERS_DEFAULT):
        print(f"{'=' * 50}")
        print(f"Solver: {solver['name']}")
        print(f"{'=' * 50}\n")
        for snr_idx, snr in enumerate(snr_values):
            base_seed = 24 + snr_idx * 1000
            lbl = f"SNR={snr:.2f} [{solver['name']}]"
            result = run_mc_trex(
                n=n, p=p, s=s, snr=snr,
                solver_desc=solver,
                base_ctrl_dict=base_ctrl,
                tFDR=tFDR,
                base_seed=base_seed,
                num_MC=num_MC,
                max_workers=num_workers,
                fixed_support=None,   # redraw per trial via make_support_random
                amplitude=1.0,
                rnd_coef=rnd_coef,
                track_L_T=True,
                label=lbl,
            )
            fdr_mat[si, snr_idx]   = result["mean_FDR"]
            tpr_mat[si, snr_idx]   = result["mean_TPR"]
            avg_L_mat[si, snr_idx] = result["mean_L"]
            avg_T_mat[si, snr_idx] = result["mean_T"]
        print()

    solver_names = [s["name"] for s in SOLVERS_DEFAULT]
    stem = f"demo_trex_03_mc_sim_variable_support_results_n{n}_p{p}_stagnation_window_{stagnant_sw}"
    save_mc_results(
        _OUT_DIR, stem, num_MC,
        solver_names, snr_values,
        fdr_mat, tpr_mat, avg_L_mat, avg_T_mat,
    )


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    trx_03_mc_sim_variable_support()
    print("\ndemo_trex_03_mc_sim_variable_support.py complete.")
