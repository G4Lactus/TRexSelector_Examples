# ==============================================================================
# demo_trex_scr_04_biobank_inmem.py
# ==============================================================================
#
# Biobank Screen-TRex Selector ("Algorithm 1"), in-memory variant.
#
# The biobank workflow screens MANY phenotypes against ONE shared design matrix
# X. For each phenotype it routes adaptively:
#   1. try Screen-TRex Ordinary; if its estimated FDR lands in the acceptable
#      window [lower_bound_FDR, upper_bound_FDR] it is accepted;
#   2. otherwise try Screen-TRex Bootstrap-CI;
#   3. otherwise fall back to the classical T-Rex selector at target_FDR_trex.
# The fraction of phenotypes routed to each method is the "Usage %".
#
# Result objects (BiobankScreenTRexResult) expose:
#   method_used ("Screen-TRex (ordinary)" / "Screen-TRex (bootstrap-CI)" /
#                "T-Rex (fallback)"),
#   estimated_FDR, estimated_FDR_screen_ordinary, estimated_FDR_screen_bootstrap,
#   selected_indices (0-based), used_fallback_trex, phenotype_index.
# A 1-D response returns one result; a m-D response returns a list of results.
#
# Mirrors:
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_04_biobank_screen_trex_inmem/
#       demo_trex_scr_04_biobank_screen_trex_inmem.cpp
#   R/trex_selector_methods/trex_screening/demo_trex_scr_04_biobank_inmem/
#       demo_trex_scr_04_biobank_inmem.R
#
# Demo content:
#   Part 1: single phenotype - basic Algorithm 1 functionality.
#   Part 2: multiple phenotypes (q=5) sharing one X, per-phenotype routing table.
#   Part 3: MC SNR sweep (single phenotype): per-method Usage %, FDR, TPR, Est.FDR.
#   Part 4: MC SNR sweep (multiple phenotypes): same per-method metrics.
#
# The C++ demo uses num_MC = 500; here we use smaller counts. Support indices
# are 0-based.
# ==============================================================================

import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np
import trex_selector_neo as ts
from trex_selector_neo.utils import compute_fdp, compute_tpp
from trex_scr_common import dgp_iid_snr, _make_response
from trex_scr_bb_common import (BIOBANK_METHODS, make_biobank_ctrl,
                                print_biobank_single, print_biobank_table,
                                print_biobank_summary, save_and_print_biobank_mc)

OUT_DIR = os.path.join(_THIS_DIR, "simulation_results", "data")

RUN_PART_1 = True
RUN_PART_2 = True
RUN_PART_3 = True
RUN_PART_4 = True

N, P = 300, 1000
NUM_MC_SINGLE = 40   # C++ uses 500; reduced for a quick runtime.
NUM_MC_MULTI = 15    # C++ uses 500; reduced (q=5 phenotypes per run).


def run_biobank(X, Y, seed=42, R_boot=1000, use_mmap=False):
    """Build a biobank selector, run it, return the select() result(s)."""
    ctrl = make_biobank_ctrl(use_mmap=use_mmap, R_boot=R_boot,
                             lower_bound_FDR=0.05, upper_bound_FDR=0.15,
                             target_FDR_trex=0.10)
    return ts.TRexBiobankScreeningSelector(X, Y, bio_ctrl=ctrl, seed=seed,
                                           verbose=False).select()


# ==============================================================================
# Part 1: single phenotype
# ==============================================================================

def part_1():
    print("\n" + "=" * 70)
    print("Part 1: Biobank Screen-TRex - Single Phenotype (in-memory)")
    print("=" * 70)
    true_support = [4, 27, 42, 149, 398]   # 0-based
    snr = 1.0
    print(f"DGP: n = {N}, p = {P}, p1 = {len(true_support)}, SNR = {snr:.1f}")
    dat = dgp_iid_snr(N, P, true_support, [1.0] * len(true_support), snr, seed=123)
    res = run_biobank(dat["X"], dat["y"])
    print_biobank_single(res, dat["true_support"])


# ==============================================================================
# Part 2: multiple phenotypes sharing one design matrix
# ==============================================================================

def part_2():
    print("\n" + "=" * 70)
    print("Part 2: Biobank Screen-TRex - Multiple Phenotypes (in-memory)")
    print("=" * 70)
    q = 5
    snr_values = [1.0, 2.0, 5.0, 10.0, 20.0]
    sup_sizes = [5, 6, 5, 7, 5]

    # Shared X (fixed seed), per-phenotype response with its own support/SNR.
    rng_x = np.random.default_rng(4212)
    X = rng_x.standard_normal((N, P))
    true_supports = []
    Y = np.zeros((N, q))
    rng_sup = np.random.default_rng(777)
    for k in range(q):
        sup = np.sort(rng_sup.choice(P, size=sup_sizes[k], replace=False))
        true_supports.append(sup)
        beta = np.zeros(P)
        beta[sup] = 1.0
        Y[:, k] = _make_response(X, beta, snr_values[k], rng_sup)
        print(f"Phenotype {k}: p1 = {sup_sizes[k]}, SNR = {snr_values[k]:.2f}")

    print(f"\nRunning Algorithm 1 for {q} phenotypes ...")
    results = run_biobank(X, Y, R_boot=500)
    print_biobank_table(results, true_supports)
    print_biobank_summary(results)


# ==============================================================================
# MC accumulation helpers
# ==============================================================================

def _new_maps(ns):
    def z():
        return {m: [0.0] * ns for m in BIOBANK_METHODS}
    return {"fdp": z(), "tpp": z(), "est": z(), "usage": z()}


def _accumulate(maps, res, true_support, si, tFDR):
    m = res.method_used
    sel = list(res.selected_indices)
    tsup = list(true_support)
    maps["fdp"][m][si] += compute_fdp(sel, tsup)
    maps["tpp"][m][si] += compute_tpp(sel, tsup)
    maps["usage"][m][si] += 1
    maps["est"]["Screen-TRex (ordinary)"][si] += res.estimated_FDR_screen_ordinary
    maps["est"]["Screen-TRex (bootstrap-CI)"][si] += res.estimated_FDR_screen_bootstrap
    maps["est"]["T-Rex (fallback)"][si] += tFDR


def _finalize_and_save(maps, snr_values, num_MC, cond_total, stem):
    """Normalize accumulators and write the table.

    cond_total is the denominator for est.FDR / usage (total phenotype screenings).
    """
    ns = len(snr_values)
    for m in BIOBANK_METHODS:
        for i in range(ns):
            if maps["usage"][m][i] > 0:
                maps["fdp"][m][i] /= maps["usage"][m][i]     # conditional
                maps["tpp"][m][i] /= maps["usage"][m][i]
            maps["est"][m][i] /= cond_total                   # unconditional
            maps["usage"][m][i] /= cond_total                 # fraction
    save_and_print_biobank_mc(num_MC, OUT_DIR, stem, snr_values, BIOBANK_METHODS,
                              maps["fdp"], maps["tpp"], maps["est"], maps["usage"])


# ==============================================================================
# Part 3: MC SNR sweep - single phenotype
# ==============================================================================

def part_3():
    print("\n" + "=" * 70)
    print("Part 3: Biobank Screen-TRex - MC SNR sweep, single phenotype (in-memory)")
    print(f"n = {N}, p = {P}, s = 10, num_MC = {NUM_MC_SINGLE}")
    print("=" * 70)
    support_size, tFDR = 10, 0.10
    snr_values = [0.1, 0.5, 1.0, 2.0, 5.0]
    maps = _new_maps(len(snr_values))

    for si, snr in enumerate(snr_values):
        for mc in range(NUM_MC_SINGLE):
            seed = 1000 * si + mc + 1
            rng = np.random.default_rng(seed + 500000)
            sup = rng.choice(P, size=support_size, replace=False)
            dat = dgp_iid_snr(N, P, sup, [1.0] * support_size, snr, seed)
            res = run_biobank(dat["X"], dat["y"], seed=-1, R_boot=500)
            _accumulate(maps, res, dat["true_support"], si, tFDR)
        print(f"  SNR {snr:.2f} -- completed {NUM_MC_SINGLE} runs.")
    _finalize_and_save(maps, snr_values, NUM_MC_SINGLE, NUM_MC_SINGLE,
                       "d04_biobank_screen_trex_mc_snr_n300_p1000_s10_inmem")


# ==============================================================================
# Part 4: MC SNR sweep - multiple phenotypes
# ==============================================================================

def part_4():
    print("\n" + "=" * 70)
    print("Part 4: Biobank Screen-TRex - MC SNR sweep, multiple phenotypes (in-memory)")
    print(f"n = {N}, p = {P}, q = 5, s = 5, num_MC = {NUM_MC_MULTI}")
    print("=" * 70)
    q, support_size, tFDR = 5, 5, 0.10
    snr_values = [0.1, 0.5, 1.0, 2.0, 5.0]
    maps = _new_maps(len(snr_values))

    for si, snr in enumerate(snr_values):
        for mc in range(NUM_MC_MULTI):
            x_seed = 5000 * si + mc + 1
            rng_x = np.random.default_rng(x_seed)
            X = rng_x.standard_normal((N, P))          # shared X per MC run
            true_supports = []
            Y = np.zeros((N, q))
            for k in range(q):
                rng_k = np.random.default_rng(
                    77 * 100000 + q * (NUM_MC_MULTI * si + mc) + k)
                sup = np.sort(rng_k.choice(P, size=support_size, replace=False))
                true_supports.append(sup)
                beta = np.zeros(P)
                beta[sup] = 1.0
                Y[:, k] = _make_response(X, beta, snr, rng_k)
            results = run_biobank(X, Y, seed=-1, R_boot=500)
            for k in range(q):
                _accumulate(maps, results[k], true_supports[k], si, tFDR)
        print(f"  SNR {snr:.2f} -- completed {NUM_MC_MULTI} runs "
              f"({q} phenotypes each).")
    total = q * NUM_MC_MULTI
    _finalize_and_save(maps, snr_values, NUM_MC_MULTI, total,
                       "d04_biobank_screen_trex_mc_multi_n300_p1000_q5_s5_inmem")


# ==============================================================================
# Main
# ==============================================================================

def main():
    if RUN_PART_1:
        part_1()
    if RUN_PART_2:
        part_2()
    if RUN_PART_3:
        part_3()
    if RUN_PART_4:
        part_4()
    print("\ndemo_trex_scr_04_biobank_inmem.py complete.")


if __name__ == "__main__":
    main()
