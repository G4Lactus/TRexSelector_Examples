# ==============================================================================
# demo_trex_da_06_bt_ht_block.py
# ==============================================================================
#
# DA-TRex (BT aggregation) against heavy-tailed Student-t(nu) block data.
# Predictors are AR(1)-block Gaussian converted to multivariate Student-t via a
# scale mixture. Two noise scenarios are run for every sweep:
#   Scenario 1 — Gaussian noise;  Scenario 2 — Student-t(nu) noise.
#
#   Part 1: SNR sweep.      Part 4: M sweep.
#   Part 2: rho sweep.      Part 5: tFDR sweep.
#   Part 3: Q sweep.        Part 6: linkage sweep.
#
# Each non-linkage sweep also loops the HAC linkage over {single, complete,
# average}. Mirrors R/.../demo_trex_da_06_bt_ht_block.R and
# cpp/.../demo_trex_da_06_mc_sim_bt_ht_block_sweeps (cpp uses a Cholesky heavy-
# tail DGP; the intent matches). All parts active — reduce num_MC to preview.
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from dgp_generators import dgp_ht_block_snr  # noqa: F401 (resolved in workers)
from da_sim_common import (make_support_bt_one_per_block, run_mc,
                           print_table, print_table_multi)

RUN_PART_1 = True
RUN_PART_2 = True
RUN_PART_3 = True
RUN_PART_4 = True
RUN_PART_5 = True
RUN_PART_6 = True

BASE = dict(n=150, M=5, Q=5, amplitude=1.0, rho=0.8, nu=3.0, snr=2.0, tFDR=0.2,
            K=20, num_MC=200, seed=2026)
LINKAGES = ("single", "complete", "average")
SCENARIOS = ((False, "Scenario 1: Gaussian Noise"), (True, "Scenario 2: Heavy Noise"))


def _mc(p, support, n_blocks, rho, snr, nu, heavy_noise, tFDR, lnk):
    return run_mc(BASE["n"], p, support, BASE["amplitude"], snr, tFDR, BASE["K"],
                  BASE["num_MC"], BASE["seed"], "dgp_ht_block_snr",
                  {"n_blocks": n_blocks, "rho": rho, "nu": nu, "heavy_noise": heavy_noise},
                  "BT", {"hc_linkage": lnk}, s=n_blocks)


def _hdr(title):
    print("=" * 70)
    print(f"  {title}  |  n={BASE['n']}  M={BASE['M']}  Q={BASE['Q']}  rho={BASE['rho']:.1f}"
          f"  nu={BASE['nu']:.0f}  {BASE['num_MC']} MC")
    print("=" * 70 + "\n")


# ==============================================================================
# Parts 1-5: sweeps, each looping scenario x linkage
# ==============================================================================

def part_1_snr_sweep():
    _hdr("HT-block — SNR sweep")
    p = BASE["M"] * BASE["Q"]
    support = make_support_bt_one_per_block(p, BASE["M"])
    for h_noise, scen in SCENARIOS:
        for lnk in LINKAGES:
            print(f'\n  [{scen}]  hc_linkage = "{lnk}"')
            rows = []
            for snr in [0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0]:
                print(f"  [SNR = {snr:.2f}] running {BASE['num_MC']} MC trials ...")
                rows.append(dict(SNR=snr, **_mc(p, support, BASE["M"], BASE["rho"],
                                                snr, BASE["nu"], h_noise, BASE["tFDR"], lnk)))
            print_table(rows, "SNR", value_fmt="{:<6.2f}", width=6)


def part_2_rho_sweep():
    _hdr("HT-block — rho sweep")
    p = BASE["M"] * BASE["Q"]
    support = make_support_bt_one_per_block(p, BASE["M"])
    for h_noise, scen in SCENARIOS:
        for lnk in LINKAGES:
            print(f'\n  [{scen}]  hc_linkage = "{lnk}"')
            rows = []
            for rho in [round(0.1 * i, 1) for i in range(10)]:
                print(f"  [rho = {rho:.1f}] running {BASE['num_MC']} MC trials ...")
                rows.append(dict(rho=rho, **_mc(p, support, BASE["M"], rho,
                                                BASE["snr"], BASE["nu"], h_noise, BASE["tFDR"], lnk)))
            print_table(rows, "rho", value_fmt="{:<6.1f}", width=6)


def part_3_q_sweep():
    _hdr("HT-block — Q sweep")
    for h_noise, scen in SCENARIOS:
        for lnk in LINKAGES:
            print(f'\n  [{scen}]  hc_linkage = "{lnk}"')
            rows = []
            for q in range(5, 51, 5):
                p = BASE["M"] * q
                support = make_support_bt_one_per_block(p, BASE["M"])
                print(f"  [Q = {q}, p = {p}] running {BASE['num_MC']} MC trials ...")
                rows.append(dict(Q=q, p=p, **_mc(p, support, BASE["M"], BASE["rho"],
                                                 BASE["snr"], BASE["nu"], h_noise, BASE["tFDR"], lnk)))
            print_table_multi(rows, ["Q", "p"], {"Q": "{:<4d}", "p": "{:<6d}"}, {"Q": 4, "p": 6})


def part_4_m_sweep():
    _hdr("HT-block — M sweep")
    for h_noise, scen in SCENARIOS:
        for lnk in LINKAGES:
            print(f'\n  [{scen}]  hc_linkage = "{lnk}"')
            rows = []
            for m in [1, 3, 5, 10, 15, 20, 30]:
                p = m * BASE["Q"]
                support = make_support_bt_one_per_block(p, m)
                print(f"  [M = {m}, p = {p}, s = {m}] running {BASE['num_MC']} MC trials ...")
                rows.append(dict(M=m, p=p, s=m, **_mc(p, support, m, BASE["rho"],
                                                      BASE["snr"], BASE["nu"], h_noise, BASE["tFDR"], lnk)))
            print_table_multi(rows, ["M", "p", "s"],
                              {"M": "{:<4d}", "p": "{:<6d}", "s": "{:<4d}"}, {"M": 4, "p": 6, "s": 4})


def part_5_tfdr_sweep():
    _hdr("HT-block — tFDR sweep")
    p = BASE["M"] * BASE["Q"]
    support = make_support_bt_one_per_block(p, BASE["M"])
    for h_noise, scen in SCENARIOS:
        for lnk in LINKAGES:
            print(f'\n  [{scen}]  hc_linkage = "{lnk}"')
            rows = []
            for tfdr in [round(0.05 * i, 2) for i in range(1, 11)]:
                print(f"  [tFDR = {tfdr:.2f}] running {BASE['num_MC']} MC trials ...")
                rows.append(dict(tFDR=tfdr, **_mc(p, support, BASE["M"], BASE["rho"],
                                                  BASE["snr"], BASE["nu"], h_noise, tfdr, lnk)))
            print_table(rows, "tFDR", value_fmt="{:<6.2f}", width=6)


# ==============================================================================
# Part 6: linkage sweep (linkage is the swept parameter), per scenario
# ==============================================================================

def part_6_linkage_sweep():
    _hdr("HT-block — linkage sweep")
    p = BASE["M"] * BASE["Q"]
    support = make_support_bt_one_per_block(p, BASE["M"])
    for h_noise, scen in SCENARIOS:
        print(f"\n  [{scen}]")
        rows = []
        for lnk in LINKAGES:
            print(f'  [linkage = {lnk}] running {BASE["num_MC"]} MC trials ...')
            rows.append(dict(Linkage=lnk, **_mc(p, support, BASE["M"], BASE["rho"],
                                                BASE["snr"], BASE["nu"], h_noise, BASE["tFDR"], lnk)))
        print_table(rows, "Linkage", value_fmt="{:<10s}", width=10)


if __name__ == "__main__":
    if RUN_PART_1:
        part_1_snr_sweep()
    if RUN_PART_2:
        part_2_rho_sweep()
    if RUN_PART_3:
        part_3_q_sweep()
    if RUN_PART_4:
        part_4_m_sweep()
    if RUN_PART_5:
        part_5_tfdr_sweep()
    if RUN_PART_6:
        part_6_linkage_sweep()
    print("\ndemo_trex_da_06_bt_ht_block.py complete.")
