# ==============================================================================
# demo_trex_da_08_groups.py
# ==============================================================================
#
# DA-TRex MC simulation for the prior-groups DGP. The predictors carry a known
# 3-level hierarchical group structure (groups of 10 / 50 / 250); the DA method
# PRIOR_GROUPS uses that structure directly. Data is drawn from a multi-level
# latent-factor model with Toeplitz leaf blocks.
#
#   Part 1: SNR sweep (n=300, p=1000, s=10, Random support).
#
# Mirrors cpp/.../demo_trex_da_08_mc_sim_groups. There is no R counterpart: the
# TRexSelectorNeo R binding's trex_da_control() exposes no prior_groups /
# rho_grid_labels parameter, so the prior-groups method cannot be driven from R
# (a binding gap, to be reported upstream). All parts active.
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from dgp_generators import dgp_groups_toeplitz_leaf  # noqa: F401 (resolved in workers)
from da_sim_common import run_mc, print_table

RUN_PART_1 = True

BASE = dict(n=300, p=1000, s=10, amplitude=1.0, tFDR=0.2, K=20, num_MC=200,
            seed=2026, max_gap=20)

# 3-level hierarchy: groups of 10, 50, 250 (fine -> coarse).
RHO_LEVELS = [0.55, 0.25, 0.10]
PHI_LEAF = 0.50


def _group_structure(p):
    return [[j // 10 for j in range(p)],
            [j // 50 for j in range(p)],
            [j // 250 for j in range(p)]]


def part_1_snr_sweep():
    p = BASE["p"]
    groups = _group_structure(p)
    print("=" * 70)
    print(f"  Prior-groups + Toeplitz leaf — SNR sweep  |  n={BASE['n']}  p={p}"
          f"  s={BASE['s']}  {BASE['num_MC']} MC")
    print(f"  3-level hierarchy (groups of 10 / 50 / 250), rho_levels={RHO_LEVELS},"
          f" phi_leaf={PHI_LEAF}")
    print("=" * 70 + "\n")

    rows = []
    for snr in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        print(f"  [SNR = {snr:.2f}] running {BASE['num_MC']} MC trials ...")
        r = run_mc(BASE["n"], p, None, BASE["amplitude"], snr, BASE["tFDR"],
                   BASE["K"], BASE["num_MC"], BASE["seed"], "dgp_groups_toeplitz_leaf",
                   {"groups": groups, "rho_levels": RHO_LEVELS, "phi": PHI_LEAF},
                   "PRIOR_GROUPS", {"prior_groups": groups, "rho_grid_labels": RHO_LEVELS},
                   s=BASE["s"])
        rows.append(dict(SNR=snr, **r))
    print("\n  Results — SNR sweep:")
    print_table(rows, "SNR", value_fmt="{:<8.2f}", width=8)


if __name__ == "__main__":
    if RUN_PART_1:
        part_1_snr_sweep()
    print("\ndemo_trex_da_08_groups.py complete.")
