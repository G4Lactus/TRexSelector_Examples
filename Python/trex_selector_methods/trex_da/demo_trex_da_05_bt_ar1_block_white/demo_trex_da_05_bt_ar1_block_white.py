# ==============================================================================
# demo_trex_da_05_bt_ar1_block_white.py
# ==============================================================================
#
# DA-TRex (BT aggregation) MC sweeps for the AR(1)-block + white-noise DGP:
# p_total predictors, of which p_ar = M*Q form M AR(1) blocks (structured) and
# p_white = p_total - p_ar are i.i.d. N(0,1) white noise. p_total is fixed at
# 1000 throughout; one active predictor per AR block. Each sweep loops the HAC
# linkage over {single, complete, average}.
#
#   Part 1: SNR sweep.
#   Part 2: rho sweep.
#   Part 3: Q sweep (p_ar = M*Q varies, p_white = p_total - p_ar).
#   Part 4: M sweep (p_ar, p_white, s = M all vary).
#
# Mirrors R/trex_selector_methods/trex_da/demo_trex_da_05_bt_ar1_block_white.R
# and cpp/.../demo_trex_da_05_mc_sim_bt_ar1_block_sweeps. All parts active (as in
# R); each is heavy — reduce num_MC to preview.
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from dgp_generators import dgp_ar1_block_white_snr  # noqa: F401 (resolved in workers)
from da_sim_common import (make_support_bt_one_per_block, run_mc,
                           print_table, print_table_multi)

RUN_PART_1 = True
RUN_PART_2 = True
RUN_PART_3 = True
RUN_PART_4 = True

BASE = dict(n=300, p_total=1000, M=5, Q=5, amplitude=1.0, rho=0.7, snr=2.0,
            tFDR=0.2, K=20, num_MC=200, seed=2026)
LINKAGES = ("single", "complete", "average")


def _mc(p_ar, p_white, support, n_blocks, rho, snr, lnk):
    return run_mc(BASE["n"], p_ar + p_white, support, BASE["amplitude"], snr,
                  BASE["tFDR"], BASE["K"], BASE["num_MC"], BASE["seed"],
                  "dgp_ar1_block_white_snr",
                  {"p_ar": p_ar, "p_white": p_white, "n_blocks": n_blocks, "rho": rho},
                  "BT", {"hc_linkage": lnk}, s=n_blocks)


def part_1_snr_sweep():
    p_ar = BASE["M"] * BASE["Q"]
    p_white = BASE["p_total"] - p_ar
    support = make_support_bt_one_per_block(p_ar, BASE["M"])
    print("=" * 70)
    print(f"  AR(1)+white — SNR sweep  |  n={BASE['n']}  p_total={BASE['p_total']}"
          f"  M={BASE['M']}  Q={BASE['Q']}  p_ar={p_ar}  p_white={p_white}  {BASE['num_MC']} MC")
    print("=" * 70 + "\n")
    for lnk in LINKAGES:
        print(f'\n  --- hc_linkage = "{lnk}" ---')
        rows = []
        for snr in [0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0]:
            print(f"  [SNR = {snr:.2f}] running {BASE['num_MC']} MC trials ...")
            rows.append(dict(SNR=snr, **_mc(p_ar, p_white, support, BASE["M"], BASE["rho"], snr, lnk)))
        print(f'\n  hc_linkage = "{lnk}"')
        print_table(rows, "SNR", value_fmt="{:<6.2f}", width=6)


def part_2_rho_sweep():
    p_ar = BASE["M"] * BASE["Q"]
    p_white = BASE["p_total"] - p_ar
    support = make_support_bt_one_per_block(p_ar, BASE["M"])
    print("=" * 70)
    print(f"  AR(1)+white — rho sweep  |  n={BASE['n']}  p_total={BASE['p_total']}"
          f"  M={BASE['M']}  Q={BASE['Q']}  SNR={BASE['snr']:.1f}  {BASE['num_MC']} MC")
    print("=" * 70 + "\n")
    for lnk in LINKAGES:
        print(f'\n  --- hc_linkage = "{lnk}" ---')
        rows = []
        for rho in [round(0.1 * i, 1) for i in range(10)]:
            print(f"  [rho = {rho:.1f}] running {BASE['num_MC']} MC trials ...")
            rows.append(dict(rho=rho, **_mc(p_ar, p_white, support, BASE["M"], rho, BASE["snr"], lnk)))
        print(f'\n  hc_linkage = "{lnk}"')
        print_table(rows, "rho", value_fmt="{:<6.1f}", width=6)


def part_3_q_sweep():
    print("=" * 70)
    print(f"  AR(1)+white — Q sweep  |  n={BASE['n']}  p_total={BASE['p_total']}"
          f"  M={BASE['M']}  rho={BASE['rho']:.1f}  SNR={BASE['snr']:.1f}  {BASE['num_MC']} MC")
    print("=" * 70 + "\n")
    for lnk in LINKAGES:
        print(f'\n  --- hc_linkage = "{lnk}" ---')
        rows = []
        for q in range(5, 51, 5):
            p_ar = BASE["M"] * q
            p_white = BASE["p_total"] - p_ar
            support = make_support_bt_one_per_block(p_ar, BASE["M"])
            print(f"  [Q = {q}, p_ar = {p_ar}] running {BASE['num_MC']} MC trials ...")
            rows.append(dict(Q=q, p_ar=p_ar, p_white=p_white,
                             **_mc(p_ar, p_white, support, BASE["M"], BASE["rho"], BASE["snr"], lnk)))
        print(f'\n  hc_linkage = "{lnk}"')
        print_table_multi(rows, ["Q", "p_ar", "p_white"],
                          {"Q": "{:<4d}", "p_ar": "{:<6d}", "p_white": "{:<8d}"},
                          {"Q": 4, "p_ar": 6, "p_white": 8})


def part_4_m_sweep():
    print("=" * 70)
    print(f"  AR(1)+white — M sweep  |  n={BASE['n']}  p_total={BASE['p_total']}"
          f"  Q={BASE['Q']}  rho={BASE['rho']:.1f}  SNR={BASE['snr']:.1f}  {BASE['num_MC']} MC")
    print("=" * 70 + "\n")
    for lnk in LINKAGES:
        print(f'\n  --- hc_linkage = "{lnk}" ---')
        rows = []
        for m in [1, 3, 5, 10, 15, 20, 30]:
            p_ar = m * BASE["Q"]
            p_white = BASE["p_total"] - p_ar
            support = make_support_bt_one_per_block(p_ar, m)
            print(f"  [M = {m}, p_ar = {p_ar}, s = {m}] running {BASE['num_MC']} MC trials ...")
            rows.append(dict(M=m, s=m, p_ar=p_ar, p_white=p_white,
                             **_mc(p_ar, p_white, support, m, BASE["rho"], BASE["snr"], lnk)))
        print(f'\n  hc_linkage = "{lnk}"')
        print_table_multi(rows, ["M", "s", "p_ar", "p_white"],
                          {"M": "{:<4d}", "s": "{:<4d}", "p_ar": "{:<6d}", "p_white": "{:<8d}"},
                          {"M": 4, "s": 4, "p_ar": 6, "p_white": 8})


if __name__ == "__main__":
    if RUN_PART_1:
        part_1_snr_sweep()
    if RUN_PART_2:
        part_2_rho_sweep()
    if RUN_PART_3:
        part_3_q_sweep()
    if RUN_PART_4:
        part_4_m_sweep()
    print("\ndemo_trex_da_05_bt_ar1_block_white.py complete.")
