# ==============================================================================
# demo_trex_scr_02_screen_mmap.py
# ==============================================================================
#
# Screen-TRex Selector, part 02: memory-mapped workflow.
#
# Statistically identical to demo 01, but the dummy matrices (and, in Part B,
# the design matrix X) are backed by memory-mapped files rather than held fully
# in RAM. This is the regime Screen-TRex is built for: p so large that a plain
# in-core T-Rex run is impractical.
#   * use_memory_mapping = True in the T-Rex control activates the internal
#     D-mmap + solver-serialization pipeline.
#   * numpy_to_memmap() writes an in-memory X to a memory-mapped binary file so
#     X itself never has to reside fully in RAM.
#
# Mirrors:
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_02_screen_trex_mmap/
#       demo_trex_scr_02_screen_trex_mmap.cpp
#   R/trex_selector_methods/trex_screening/demo_trex_scr_02_screen_mmap/
#       demo_trex_scr_02_screen_mmap.R
#
# Demo content:
#   Part A: in-memory X + use_memory_mapping = True (D-mmap only),
#           Ordinary and Bootstrap-CI single runs.
#   Part B: fully memory-mapped X (written to a temp mmap file) + D-mmap,
#           single Ordinary run.
#   Part C: Monte Carlo SNR sweep (Ordinary vs. Bootstrap-CI) with mmap enabled.
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
from trex_selector_neo.utils import numpy_to_memmap
from trex_scr_common import (dgp_iid_snr, make_scr_trex_ctrl, make_screen_ctrl,
                             default_scr_methods, print_scr_result,
                             run_mc_screen, save_and_print_scr_mc)

OUT_DIR = os.path.join(_THIS_DIR, "simulation_results", "data")

RUN_PART_A = True
RUN_PART_B = True
RUN_PART_C = True

TRUE_SUPPORT = [27, 149, 398, 4, 42]   # 0-based, matches the C++ demo
NUM_MC = 40   # C++ uses 500; reduced for a quick runtime.


# ==============================================================================
# Part A: in-memory X + use_memory_mapping = True
# ==============================================================================

def part_a():
    n, p, snr = 300, 1000, 5.0
    for bootstrap in (False, True):
        label = "Bootstrap-CI" if bootstrap else "Ordinary"
        print("\n" + "=" * 70)
        print(f"Part A: in-memory X + D-mmap  |  {label} screening")
        print(f"n = {n}, p = {p}, s = {len(TRUE_SUPPORT)}, SNR = {snr:.1f}\n")

        dat = dgp_iid_snr(n, p, TRUE_SUPPORT, [1.0] * len(TRUE_SUPPORT), snr, seed=58)
        sc = make_screen_ctrl(trex_method="TREX", bootstrap=bootstrap)
        sel = ts.TRexScreeningSelector(
            dat["X"], dat["y"], screen_control=sc,
            trex_control=make_scr_trex_ctrl(use_memory_mapping=True),
            seed=42, verbose=False)
        res = sel.select()
        print_scr_result(sel, res, dat["true_support"])


# ==============================================================================
# Part B: fully memory-mapped X + use_memory_mapping = True
# ==============================================================================

def part_b():
    n, p, snr = 300, 1000, 5.0
    print("\n" + "=" * 70)
    print("Part B: memory-mapped X (temp mmap file) + D-mmap  |  Ordinary screening")
    print(f"n = {n}, p = {p}, s = {len(TRUE_SUPPORT)}, SNR = {snr:.1f}\n")

    dat = dgp_iid_snr(n, p, TRUE_SUPPORT, [1.0] * len(TRUE_SUPPORT), snr, seed=58)
    print("Converting X to a memory-mapped file and running Screen-TRex ...")
    fd, mmap_path = tempfile.mkstemp(suffix=".dat")
    os.close(fd)
    try:
        X_mmap = numpy_to_memmap(mmap_path, dat["X"])
        sc = make_screen_ctrl(trex_method="TREX")
        sel = ts.TRexScreeningSelector(
            X_mmap.to_numpy(), dat["y"], screen_control=sc,
            trex_control=make_scr_trex_ctrl(use_memory_mapping=True),
            seed=42, verbose=False)
        res = sel.select()
        print_scr_result(sel, res, dat["true_support"])
    finally:
        try:
            os.unlink(mmap_path)
        except OSError:
            pass


# ==============================================================================
# Part C: Monte Carlo SNR sweep with mmap enabled
# ==============================================================================

def part_c():
    print("\n" + "=" * 70)
    print("Part C: Monte Carlo SNR sweep (Ordinary vs. Bootstrap-CI), mmap enabled")
    print(f"n = 300, p = 1000, s = 10, num_MC = {NUM_MC}\n")

    n, p, support_size = 300, 1000, 10
    snr_values = [0.1, 0.5, 1.0, 2.0, 5.0]
    methods = default_scr_methods()
    method_names = [m["name"] for m in methods]

    def init():
        return {nm: [0.0] * len(snr_values) for nm in method_names}
    fdr_map, tpr_map, est_fdr_map = init(), init(), init()

    trex_ctrl = make_scr_trex_ctrl(use_memory_mapping=True)

    for m in methods:
        print("=" * 70)
        print(f"Method: {m['name']}")
        print("=" * 70)
        for si, snr in enumerate(snr_values):
            def make_data(seed, snr=snr):
                rng = np.random.default_rng(seed + 500000)
                sup = rng.choice(p, size=support_size, replace=False)
                return dgp_iid_snr(n, p, sup, [1.0] * support_size, snr, seed)
            res = run_mc_screen(
                NUM_MC, f"{m['name']}  SNR={snr:.2f}", make_data, trex_ctrl,
                screen_args={"trex_method": m["trex_method"],
                             "bootstrap": m["bootstrap"]},
                base_seed=24 + si * 1000)
            fdr_map[m["name"]][si] = res["avg_fdr"]
            tpr_map[m["name"]][si] = res["avg_tpr"]
            est_fdr_map[m["name"]][si] = res["avg_est_fdr"]

    save_and_print_scr_mc(NUM_MC, OUT_DIR, "d02_screen_trex_mmap_mc_n300_p1000",
                          snr_values, method_names, fdr_map, tpr_map, est_fdr_map)


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
    print("\ndemo_trex_scr_02_screen_mmap.py complete.")


if __name__ == "__main__":
    main()
