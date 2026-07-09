# ==============================================================================
# demo_trex_scr_06_screen_solvers.py
# ==============================================================================
#
# Screen-TRex Selector, part 06: solver-backend comparison.
#
# Screening thresholds the dummy-based voting proportion Phi_j, but Phi itself
# is produced by an underlying T-Rex solver. This demo asks whether the choice
# of solver backend changes screening performance:
#   * TLARS - T-Rex Least-Angle Regression Solver (the default),
#   * TAFS  - T-Rex Adaptive Forward Selection (rho_afs = 0.3),
#   * TOMP  - T-Rex Orthogonal Matching Pursuit,
# each combined with Ordinary and Bootstrap-CI screening.
#
# Mirrors:
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_trex_solvers/
#       demo_trex_scr_06_screen_trex_solvers.cpp
#   R/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_solvers/
#       demo_trex_scr_06_screen_solvers.R
#
# Demo content:
#   Part 1: MC SNR sweep - three solvers, Ordinary screening only (3 series).
#   Part 2: MC SNR sweep - three solvers x {Ordinary, Bootstrap-CI} (6 series).
#
# The C++ demo uses num_MC = 500; here we use a smaller count. Support indices
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
from trex_scr_common import (dgp_iid_snr, make_scr_trex_ctrl, run_mc_screen,
                             save_and_print_scr_mc)

OUT_DIR = os.path.join(_THIS_DIR, "simulation_results")

RUN_PART_1 = True
RUN_PART_2 = True

N, P, SUPPORT_SIZE = 300, 1000, 10
SNR_GRID = [0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0]
NUM_MC = 20   # C++ uses 500; reduced for a quick runtime.


def sv(name, solver, rho_afs=0.0, bootstrap=False):
    """A solver x screening-variant descriptor."""
    return {"name": name, "solver": solver, "rho_afs": rho_afs, "bootstrap": bootstrap}


def run_solver_sweep(title, variants, file_stem):
    """MC sweep over a list of solver x method variants (SNR grid)."""
    print("\n" + "=" * 70)
    print(title)
    print(f"High-dimensional: n = {N}, p = {P}, s = {SUPPORT_SIZE}, num_MC = {NUM_MC}")
    print("=" * 70 + "\n")

    method_names = [v["name"] for v in variants]

    def init():
        return {nm: [0.0] * len(SNR_GRID) for nm in method_names}
    fdr_map, tpr_map, est_fdr_map = init(), init(), init()

    for v in variants:
        print("=" * 70)
        print(f"Variant: {v['name']}")
        print("=" * 70)
        # TAFS/TOMP: enable stagnation stop, as in the C++ make_trex_ctrl_for().
        extra = {"tloop_stagnation_stop": True, "tloop_max_stagnant_steps": 5}
        if v["rho_afs"]:
            extra["rho_afs"] = v["rho_afs"]
        trex_ctrl = make_scr_trex_ctrl(solver=v["solver"], **extra)
        for si, snr in enumerate(SNR_GRID):
            def make_data(seed, snr=snr):
                rng = np.random.default_rng(seed + 500000)
                sup = rng.choice(P, size=SUPPORT_SIZE, replace=False)
                return dgp_iid_snr(N, P, sup, [1.0] * SUPPORT_SIZE, snr, seed)
            res = run_mc_screen(
                NUM_MC, f"{v['name']}  SNR={snr:.2f}", make_data, trex_ctrl,
                screen_args={"trex_method": "TREX", "bootstrap": v["bootstrap"]},
                base_seed=24 + si * 1000)
            fdr_map[v["name"]][si] = res["avg_fdr"]
            tpr_map[v["name"]][si] = res["avg_tpr"]
            est_fdr_map[v["name"]][si] = res["avg_est_fdr"]

    save_and_print_scr_mc(NUM_MC, OUT_DIR, file_stem, SNR_GRID,
                          method_names, fdr_map, tpr_map, est_fdr_map)


# ==============================================================================
# Part 1: three solvers, Ordinary screening only
# ==============================================================================

def part_1():
    variants = [
        sv("TLARS (Ordinary)",    "TLARS"),
        sv("TAFS-0.3 (Ordinary)", "TAFS", rho_afs=0.3),
        sv("TOMP (Ordinary)",     "TOMP")]
    run_solver_sweep("Part 1: Solver comparison - Ordinary Screen-TRex (3 series)",
                     variants, "d06_screen_trex_solvers_mc_snr_n300_p1000_s10")


# ==============================================================================
# Part 2: three solvers x {Ordinary, Bootstrap-CI}
# ==============================================================================

def part_2():
    variants = [
        sv("TLARS (Ordinary)",     "TLARS"),
        sv("TAFS-0.3 (Ordinary)",  "TAFS", rho_afs=0.3),
        sv("TOMP (Ordinary)",      "TOMP"),
        sv("TLARS (Bootstrap)",    "TLARS", bootstrap=True),
        sv("TAFS-0.3 (Bootstrap)", "TAFS", rho_afs=0.3, bootstrap=True),
        sv("TOMP (Bootstrap)",     "TOMP", bootstrap=True)]
    run_solver_sweep("Part 2: Solver x method comparison (6 series)",
                     variants, "d06_screen_trex_solver_method_mc_snr_n300_p1000_s10")


# ==============================================================================
# Main
# ==============================================================================

def main():
    if RUN_PART_1:
        part_1()
    if RUN_PART_2:
        part_2()
    print("\ndemo_trex_scr_06_screen_solvers.py complete.")


if __name__ == "__main__":
    main()
