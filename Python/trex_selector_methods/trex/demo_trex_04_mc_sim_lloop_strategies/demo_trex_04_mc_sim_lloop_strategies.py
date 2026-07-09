# ==============================================================================
# demo_trex_04_mc_sim_lloop_strategies.py
# ==============================================================================
#
# Classical T-Rex selector demo — L-loop strategy comparison.
# Mirrors R/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies.R and
# cpp/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies.cpp.
#
# Demo content:
#
#    STANDARD          — Fresh i.i.d. dummy matrix at each L-loop iteration
#                        (conservative gold standard; highest memory cost).
#    HCONCAT           — Horizontally expand dummy columns from the same base draw.
#    PERMUTATION       — Re-use the base dummy matrix via permutations.
#    PERMUTATION_DIRECT— Seed-based permutations; no base matrix in memory.
#    DIRECT            — Seed-based i.i.d. draws; no base matrix in memory.
#    SKIPL             — Skip the L-loop entirely; always uses
#                        L = max_dummy_multiplier (no adaptive L calibration).
#
#  Scenario:
#  High-dimensional (n = 300, p = 1000, s = 10).
#  SNR sweep: {0.1, 0.2, ..., 2.0, 5.0}  (21 values).
#  Two support scenarios:
#    random support — redrawn per trial;
#    block support  — contiguous block {0, 1, ..., s-1}  (0-based).
#  Reports averaged FDR / TPR / Avg L / Avg T per L-strategy x SNR.
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

from trex_sim_common import run_mc_trex, save_mc_results

# ==============================================================================
# Global parameters
# ==============================================================================

_OUT_DIR     = os.path.join(_THIS_DIR, "simulation_results")
_NUM_WORKERS = 6
_NUM_MC      = 25

# ==============================================================================
# L-loop strategies to compare
# ==============================================================================
#
# name           : row label in output tables and CSV.
# strategy       : string passed to lloop_strategy in _make_trex_ctrl_from_dict.
# max_dummy_mult : only meaningful for SKIPL variants (fixed L = max_dummy_mult);
#                  adaptive strategies default to 10.

def _strat(name, strategy, max_dummy_mult=10):
    return dict(name=name, strategy=strategy, max_dummy_mult=max_dummy_mult)

L_STRATEGIES = [
    # _strat("STANDARD",           "STANDARD"),
    # _strat("HCONCAT",            "HCONCAT"),
    # _strat("PERMUTATION",        "PERMUTATION"),
    # _strat("PERMUTATION_DIRECT", "PERMUTATION_DIRECT"),
    # _strat("DIRECT",             "DIRECT"),
    # _strat("SKIPL_5p",  "SKIPL",  5),
    # _strat("SKIPL_10p", "SKIPL", 10),
    _strat("SKIPL_20p", "SKIPL", 20),
    _strat("SKIPL_50p", "SKIPL", 50),
]

# TLARS is the base solver for all L-strategy comparisons
_TLARS_SOLVER = dict(
    solver_name   = "TLARS",
    name          = "TLARS",
    lambda2       = 0.1,
    rho_afs       = 0.3,
    ncgmp_variant = 0,
)


# ==============================================================================
# MC simulation — L-loop strategy comparison
# ==============================================================================

def trx_04_lloop_strategies(num_MC=_NUM_MC, num_workers=_NUM_WORKERS,
                             rnd_coef=False, block_support=False):

    n = 300
    p = 1000
    s = 10
    tFDR       = 0.1
    snr_values = list(np.round(np.arange(0.1, 2.1, 0.1), 1)) + [5.0]  # 21 values

    support_label = "block" if block_support else "random"
    print(f"\n{'=' * 70}")
    print("MC Simulation \u2014 L-Strategy Comparison  |  TLARS")
    print("  High-dimensional (p > n)")
    print(f"  n = {n},  p = {p},  s = {s},  num_MC = {num_MC}  [{support_label} support]\n")

    # Block support: contiguous 0-based {0, 1, ..., s-1}  (fixed across all trials)
    block_supp = np.arange(s, dtype=np.intp) if block_support else None

    n_strat = len(L_STRATEGIES)
    n_snr   = len(snr_values)
    fdr_mat   = np.zeros((n_strat, n_snr))
    tpr_mat   = np.zeros((n_strat, n_snr))
    avg_L_mat = np.zeros((n_strat, n_snr))
    avg_T_mat = np.zeros((n_strat, n_snr))

    # Simulation sweep: iterate over L-loop strategies and SNR values
    for si, strat in enumerate(L_STRATEGIES):
        print(f"{'=' * 50}")
        print(f"L-strategy: {strat['name']}")
        print(f"{'=' * 50}\n")

        # Per-strategy ctrl — max_dummy_multiplier and lloop_strategy vary
        base_ctrl = dict(
            K                        = 20,
            max_dummy_multiplier     = strat["max_dummy_mult"],
            use_max_T_stop           = True,
            lloop_strategy           = strat["strategy"],
            parallel_rnd_experiments = False,   # mirrors R use_openmp=FALSE
        )

        for snr_idx, snr in enumerate(snr_values):
            base_seed = 24 + snr_idx * 1000
            lbl = f"SNR={snr:.2f} [{strat['name']}]"
            result = run_mc_trex(
                n=n, p=p, s=s, snr=snr,
                solver_desc=_TLARS_SOLVER,
                base_ctrl_dict=base_ctrl,
                tFDR=tFDR,
                base_seed=base_seed,
                num_MC=num_MC,
                max_workers=num_workers,
                fixed_support=block_supp,   # None → random per trial
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

    strat_names = [st["name"] for st in L_STRATEGIES]
    stem = (f"demo_trex_04_lloop_strategies_results_n{n}_p{p}"
            f"_{support_label}_support")
    save_mc_results(
        _OUT_DIR, stem, num_MC,
        strat_names, snr_values,
        fdr_mat, tpr_mat, avg_L_mat, avg_T_mat,
    )
    print()


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    trx_04_lloop_strategies(block_support=False)
    # trx_04_lloop_strategies(block_support=True)
    print("\ndemo_trex_04_mc_sim_lloop_strategies.py complete.")
