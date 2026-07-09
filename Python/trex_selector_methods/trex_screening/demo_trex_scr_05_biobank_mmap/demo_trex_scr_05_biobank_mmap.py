# ==============================================================================
# demo_trex_scr_05_biobank_mmap.py
# ==============================================================================
#
# Biobank Screen-TRex Selector ("Algorithm 1"), memory-mapped variant.
#
# Statistically identical to demo 04, but the design matrix X and the internal
# dummy matrices are backed by memory-mapped files (use_memory_mapping = True),
# so a biobank-scale X never has to reside fully in RAM. Everything else about
# Algorithm 1 - the Ordinary -> Bootstrap-CI -> T-Rex-fallback routing and its
# per-method Usage % - is unchanged. See demo 04 for the result-object accessors
# and method-used labels (identical here).
#
# Mirrors:
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_05_biobank_screen_trex_mmap/
#       demo_trex_scr_05_biobank_screen_trex_mmap.cpp
#   R/trex_selector_methods/trex_screening/demo_trex_scr_05_biobank_mmap/
#       demo_trex_scr_05_biobank_mmap.R
#
# Demo content:
#   Part A: single phenotype, in-memory X + use_memory_mapping = True (D-mmap).
#   Part B: single + multiple phenotypes, fully memory-mapped X (temp mmap file).
#   Part C: MC SNR sweep (single phenotype) with the mmap pipeline enabled.
#
# The C++ demo uses num_MC = 500; here we use a smaller count. Support indices
# are 0-based.
# ==============================================================================

import os
import sys
import tempfile

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np
import trex_selector_neo as ts
from trex_selector_neo.utils import compute_fdp, compute_tpp, numpy_to_memmap
from trex_scr_common import dgp_iid_snr, _make_response
from trex_scr_bb_common import (BIOBANK_METHODS, make_biobank_ctrl,
                                print_biobank_single, print_biobank_table,
                                print_biobank_summary, save_and_print_biobank_mc)

OUT_DIR = os.path.join(_THIS_DIR, "simulation_results")

RUN_PART_A = True
RUN_PART_B = True
RUN_PART_C = True

N, P = 300, 1000
NUM_MC = 30   # C++ uses 500; reduced for a quick runtime.


def run_biobank_mmap(X, Y, seed=42, R_boot=1000):
    """Run a memory-mapped biobank selector on X, Y."""
    ctrl = make_biobank_ctrl(use_mmap=True, R_boot=R_boot,
                             lower_bound_FDR=0.05, upper_bound_FDR=0.15,
                             target_FDR_trex=0.10)
    return ts.TRexBiobankScreeningSelector(X, Y, bio_ctrl=ctrl, seed=seed,
                                           verbose=False).select()


def with_temp_mmap(X_matrix, fn):
    """RAII-like temp-mmap lifecycle helper (see demo 02): write X to a temp
    mmap file, yield the mmap view, and unlink on exit even under error."""
    fd, x_path = tempfile.mkstemp(suffix=".dat")
    os.close(fd)
    try:
        X_mmap = numpy_to_memmap(x_path, X_matrix)
        return fn(X_mmap.to_numpy())
    finally:
        try:
            os.unlink(x_path)
        except OSError:
            pass


# ==============================================================================
# Part A: single phenotype, in-memory X + D-mmap
# ==============================================================================

def part_a():
    print("\n" + "=" * 70)
    print("Part A: Biobank Screen-TRex - single phenotype, in-memory X + D-mmap")
    print("=" * 70)
    true_support = [4, 27, 42, 149, 398]   # 0-based
    dat = dgp_iid_snr(N, P, true_support, [1.0] * len(true_support), snr=1.0, seed=123)
    res = run_biobank_mmap(dat["X"], dat["y"])
    print_biobank_single(res, dat["true_support"])


# ==============================================================================
# Part B: fully memory-mapped X (temp mmap file)
# ==============================================================================

def part_b():
    print("\n" + "=" * 70)
    print("Part B: Biobank Screen-TRex - memory-mapped X (single + multi phenotype)")
    print("=" * 70)

    # -- single phenotype on mmap X --
    true_support = [4, 27, 42, 149, 398]
    dat = dgp_iid_snr(N, P, true_support, [1.0] * len(true_support), snr=1.0, seed=123)
    print("Single phenotype on memory-mapped X ...")
    with_temp_mmap(dat["X"], lambda Xm:
                   print_biobank_single(run_biobank_mmap(Xm, dat["y"]),
                                        dat["true_support"]))

    # -- multiple phenotypes sharing one mmap X --
    q = 5
    snr_values = [1.0, 2.0, 5.0, 10.0, 20.0]
    sup_sizes = [5, 6, 5, 7, 5]
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

    print(f"\n{q} phenotypes on one shared memory-mapped X ...")

    def _run_multi(Xm):
        results = run_biobank_mmap(Xm, Y, R_boot=500)
        print_biobank_table(results, true_supports)
        print_biobank_summary(results)
    with_temp_mmap(X, _run_multi)


# ==============================================================================
# Part C: MC SNR sweep - single phenotype, mmap enabled
# ==============================================================================

def part_c():
    print("\n" + "=" * 70)
    print("Part C: Biobank Screen-TRex - MC SNR sweep (single phenotype), mmap enabled")
    print(f"n = {N}, p = {P}, s = 10, num_MC = {NUM_MC}")
    print("=" * 70)
    support_size, tFDR = 10, 0.10
    snr_values = [0.1, 0.5, 1.0, 2.0, 5.0]
    ns = len(snr_values)

    def z():
        return {m: [0.0] * ns for m in BIOBANK_METHODS}
    fdp, tpp, est, usage = z(), z(), z(), z()

    for si, snr in enumerate(snr_values):
        for mc in range(NUM_MC):
            seed = 1000 * si + mc + 1
            rng = np.random.default_rng(seed + 500000)
            sup = rng.choice(P, size=support_size, replace=False)
            dat = dgp_iid_snr(N, P, sup, [1.0] * support_size, snr, seed)
            res = with_temp_mmap(dat["X"], lambda Xm:
                                 run_biobank_mmap(Xm, dat["y"], seed=-1, R_boot=500))
            m = res.method_used
            sel = list(res.selected_indices)
            tsup = dat["true_support"].tolist()
            fdp[m][si] += compute_fdp(sel, tsup)
            tpp[m][si] += compute_tpp(sel, tsup)
            usage[m][si] += 1
            est["Screen-TRex (ordinary)"][si] += res.estimated_FDR_screen_ordinary
            est["Screen-TRex (bootstrap-CI)"][si] += res.estimated_FDR_screen_bootstrap
            est["T-Rex (fallback)"][si] += tFDR
        print(f"  SNR {snr:.2f} -- completed {NUM_MC} runs.")

    for m in BIOBANK_METHODS:
        for i in range(ns):
            if usage[m][i] > 0:
                fdp[m][i] /= usage[m][i]
                tpp[m][i] /= usage[m][i]
            est[m][i] /= NUM_MC
            usage[m][i] /= NUM_MC
    save_and_print_biobank_mc(NUM_MC, OUT_DIR,
                              "d05_biobank_screen_trex_mc_snr_n300_p1000_s10_mmap",
                              snr_values, BIOBANK_METHODS, fdp, tpp, est, usage)


# ==============================================================================
# Main
# ==============================================================================

def main():
    if RUN_PART_A:
        part_a()
    if RUN_PART_B:
        part_b()
    if RUN_PART_C:
        part_c()
    print("\ndemo_trex_scr_05_biobank_mmap.py complete.")


if __name__ == "__main__":
    main()
