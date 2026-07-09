# ==============================================================================
# demo_trex_da_04_bt_ar1_block.py
# ==============================================================================
#
# DA-TRex (BT aggregation) MC sweeps for the block-diagonal AR(1) DGP: p
# predictors split into M blocks of size Q, AR(1) within each block, one active
# predictor per block. Each sweep loops the HAC linkage over {single, complete,
# average}.
#
#   Part 1: SNR sweep.
#   Part 2: rho sweep.
#   Part 3: Q sweep (block size; p = M*Q varies).
#   Part 4: M sweep (number of blocks; p = M*Q and s = M vary).
#
# Mirrors R/trex_selector_methods/trex_da/demo_trex_da_04_bt_ar1_block.R and
# cpp/.../demo_trex_da_04_mc_sim_bt_ar1_block. All parts active (as in R); each
# is heavy (num_MC x sweep points x 3 linkages) — reduce num_MC to preview.
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from dgp_generators import dgp_ar1_block_snr  # noqa: F401 (resolved by name in workers)
from da_sim_common import (make_support_bt_one_per_block, run_mc,
                           print_table, print_table_multi)

RUN_PART_1 = True
RUN_PART_2 = True
RUN_PART_3 = True
RUN_PART_4 = True

BASE = dict(n=150, M=5, Q=5, amplitude=1.0, rho=0.7, snr=2.0, tFDR=0.2,
            K=20, num_MC=200, seed=2026)
LINKAGES = ("single", "complete", "average")


def _mc(p, support, n_blocks, rho, snr, lnk):
    return run_mc(BASE["n"], p, support, BASE["amplitude"], snr, BASE["tFDR"],
                  BASE["K"], BASE["num_MC"], BASE["seed"], "dgp_ar1_block_snr",
                  {"n_blocks": n_blocks, "rho": rho}, "BT", {"hc_linkage": lnk},
                  s=n_blocks)


# ==============================================================================
# Part 1: SNR sweep
# ==============================================================================

def part_1_snr_sweep():
    p = BASE["M"] * BASE["Q"]
    support = make_support_bt_one_per_block(p, BASE["M"])
    snr_grid = [0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0]
    print("=" * 70)
    print(f"  AR(1)-block — SNR sweep  |  n={BASE['n']}  M={BASE['M']}  Q={BASE['Q']}"
          f"  rho={BASE['rho']:.1f}  {BASE['num_MC']} MC")
    print("=" * 70 + "\n")
    for lnk in LINKAGES:
        print(f'\n  --- hc_linkage = "{lnk}" ---')
        rows = []
        for snr in snr_grid:
            print(f"  [SNR = {snr:.2f}] running {BASE['num_MC']} MC trials ...")
            rows.append(dict(SNR=snr, **_mc(p, support, BASE["M"], BASE["rho"], snr, lnk)))
        print(f'\n  hc_linkage = "{lnk}"')
        print_table(rows, "SNR", value_fmt="{:<8.2f}", width=8)


# ==============================================================================
# Part 2: rho sweep
# ==============================================================================

def part_2_rho_sweep():
    p = BASE["M"] * BASE["Q"]
    support = make_support_bt_one_per_block(p, BASE["M"])
    rho_grid = [round(0.1 * i, 1) for i in range(10)]
    print("=" * 70)
    print(f"  AR(1)-block — rho sweep  |  n={BASE['n']}  M={BASE['M']}  Q={BASE['Q']}"
          f"  SNR={BASE['snr']:.1f}  {BASE['num_MC']} MC")
    print("=" * 70 + "\n")
    for lnk in LINKAGES:
        print(f'\n  --- hc_linkage = "{lnk}" ---')
        rows = []
        for rho in rho_grid:
            print(f"  [rho = {rho:.1f}] running {BASE['num_MC']} MC trials ...")
            rows.append(dict(rho=rho, **_mc(p, support, BASE["M"], rho, BASE["snr"], lnk)))
        print(f'\n  hc_linkage = "{lnk}"')
        print_table(rows, "rho", value_fmt="{:<6.1f}", width=6)


# ==============================================================================
# Part 3: Q sweep (block size; p = M*Q varies)
# ==============================================================================

def part_3_q_sweep():
    q_grid = list(range(5, 51, 5))
    print("=" * 70)
    print(f"  AR(1)-block — Q sweep  |  n={BASE['n']}  M={BASE['M']}  s={BASE['M']}"
          f"  rho={BASE['rho']:.1f}  SNR={BASE['snr']:.1f}  {BASE['num_MC']} MC")
    print("=" * 70 + "\n")
    for lnk in LINKAGES:
        print(f'\n  --- hc_linkage = "{lnk}" ---')
        rows = []
        for q in q_grid:
            p = BASE["M"] * q
            support = make_support_bt_one_per_block(p, BASE["M"])
            print(f"  [Q = {q}, p = {p}] running {BASE['num_MC']} MC trials ...")
            rows.append(dict(Q=q, p=p, **_mc(p, support, BASE["M"], BASE["rho"], BASE["snr"], lnk)))
        print(f'\n  hc_linkage = "{lnk}"')
        print_table_multi(rows, ["Q", "p"], {"Q": "{:<4d}", "p": "{:<6d}"},
                          {"Q": 4, "p": 6})


# ==============================================================================
# Part 4: M sweep (number of blocks; p = M*Q and s = M vary)
# ==============================================================================

def part_4_m_sweep():
    m_grid = [1, 3, 5, 10, 15, 20, 30]
    print("=" * 70)
    print(f"  AR(1)-block — M sweep  |  n={BASE['n']}  Q={BASE['Q']}  rho={BASE['rho']:.1f}"
          f"  SNR={BASE['snr']:.1f}  {BASE['num_MC']} MC")
    print("=" * 70 + "\n")
    for lnk in LINKAGES:
        print(f'\n  --- hc_linkage = "{lnk}" ---')
        rows = []
        for m in m_grid:
            p = m * BASE["Q"]
            support = make_support_bt_one_per_block(p, m)
            print(f"  [M = {m}, p = {p}, s = {m}] running {BASE['num_MC']} MC trials ...")
            rows.append(dict(M=m, p=p, s=m, **_mc(p, support, m, BASE["rho"], BASE["snr"], lnk)))
        print(f'\n  hc_linkage = "{lnk}"')
        print_table_multi(rows, ["M", "p", "s"],
                          {"M": "{:<4d}", "p": "{:<6d}", "s": "{:<4d}"},
                          {"M": 4, "p": 6, "s": 4})


if __name__ == "__main__":
    if RUN_PART_1:
        part_1_snr_sweep()
    if RUN_PART_2:
        part_2_rho_sweep()
    if RUN_PART_3:
        part_3_q_sweep()
    if RUN_PART_4:
        part_4_m_sweep()
    print("\ndemo_trex_da_04_bt_ar1_block.py complete.")
