# ==============================================================================
# demo_ts_04_tenet.py
# ==============================================================================
#
# Demonstration of the T-ENET (Terminating Elastic Net) solver.
#
# Mirrors cpp/tsolvers/demo_ts_04_tenet/demo_ts_04_tenet.cpp.
#
# Parts:
#   Demo 1: Basic T-ENET with early stopping (internal normalization),
#           high-dimensional (p > n) and low-dimensional (n > p) regimes.
#   Demo 2: External L2 normalization (LpNormScaler) + centered y,
#           solver run with normalize=False, intercept=False.
#   Demo 3: Serialization & warm start — save a partial path, reload,
#           continue, and compare against a reference run.
#   Demo 4: Memory-mapped X/D/y via numpy_to_memmap.
#           (cpp generates p = 500,000 on disk; scaled to p = 5,000 here.)
#
# ==============================================================================

import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

from trex_selector_neo.ml_methods import LpNormScaler, NormType
from trex_selector_neo.tsolvers.lars_based import TENET_Solver
from trex_selector_neo.utils import get_max_threads, numpy_to_memmap, set_num_threads

from ts_demo_utils import (
    gen_synthetic_data,
    print_demo_config,
    print_section_header,
    print_selection,
    print_selection_quality,
)


# ==============================================================================
# Demo 1: Basic T-Elastic Net with Early Stopping
# ==============================================================================

def demo_early_stopping(high_dim, rnd_coef, T_stop):
    print_section_header("Demo 1: Basic T-Elastic Net with Early Stopping")

    n = 1000 if high_dim else 5000
    p = 5000 if high_dim else 1000
    num_dummies = 10 * p

    true_support = [27, 149, 398, 420, 4]
    true_coefs = [-0.4, -0.2, -0.8, 1.1, 2.5] if rnd_coef else [1.0] * 5
    snr = 1.0
    lambda2 = 0.5

    print("High-dimensional (p > n)" if high_dim else "Low-dimensional (n > p)")
    print_demo_config(n, p, num_dummies, T_stop, true_support, true_coefs, snr)

    X, D, y = gen_synthetic_data(n, p, num_dummies, true_support, true_coefs,
                                 snr, seed=42)

    solver = TENET_Solver(X, D, y, lambda2, normalize=True, intercept=True,
                          verbose=True)
    solver.executeStep(T_stop, early_stop=True)

    print_selection(solver, true_support)
    print_selection_quality(solver, true_support)
    print("\nT-Elastic Net Diagnostics:")
    print(f"  Removals: {solver.getNumRemovals()}")
    print(f"  Cycling ratio: {solver.getCyclingRatio():.4f}")
    print("\n")


# ==============================================================================
# Demo 2: External Normalization (LpNormScaler)
# ==============================================================================

def demo_with_external_normalizer(high_dim, rnd_coef, T_stop):
    print_section_header("Demo 2: T-Elastic Net with External Normalization")

    n = 1000 if high_dim else 5000
    p = 5000 if high_dim else 1000
    num_dummies = 10 * p

    true_support = [4, 27, 149, 398, 420]
    true_coefs = [2.5, -0.4, -0.2, -0.8, 1.1] if rnd_coef else [1.0] * 5
    snr = 1.0
    lambda2 = 0.5

    print("High-dimensional (p > n)" if high_dim else "Low-dimensional (n > p)")
    print_demo_config(n, p, num_dummies, T_stop, true_support, true_coefs, snr)

    X, D, y = gen_synthetic_data(n, p, num_dummies, true_support, true_coefs,
                                 snr, seed=42)

    # External normalization: center + unit-L2 columns for X and D, center y.
    print("Applying external L2 normalization...")
    l2scaler_X = LpNormScaler(NormType.L2, center=True)
    l2scaler_X.fit(X)
    l2scaler_X.transform_inplace(X)

    l2scaler_D = LpNormScaler(NormType.L2, center=True)
    l2scaler_D.fit(D)
    l2scaler_D.transform_inplace(D)

    y -= y.mean()

    solver = TENET_Solver(X, D, y, lambda2, normalize=False, intercept=False,
                          verbose=True)
    solver.executeStep(T_stop, early_stop=True)

    print_selection(solver, true_support)
    print_selection_quality(solver, true_support)
    print("\nT-Elastic Net Diagnostics:")
    print(f"  Removals: {solver.getNumRemovals()}")
    print(f"  Cycling ratio: {solver.getCyclingRatio():.4f}")
    print("\n")


# ==============================================================================
# Demo 3: Serialization & Warm Start
# ==============================================================================

def demo_serialization():
    print_section_header("Demo 3: T-Elastic Net Serialization & Path Consistency")

    n, p = 100, 50
    num_dummies = p
    T_stop_final = 15
    T_stop_partial = 7
    snr = 1.0

    true_support = [10, 25, 40]
    true_coefs = [2.5, -1.8, 3.2]
    lambda2 = 0.5

    print_demo_config(n, p, num_dummies, T_stop_final, true_support,
                      true_coefs, snr)

    X, D, y = gen_synthetic_data(n, p, num_dummies, true_support, true_coefs,
                                 snr, seed=42)

    # Reference solver: run until T_stop_final.
    solver_ref = TENET_Solver(X, D, y, lambda2, normalize=True, intercept=True,
                              verbose=False)
    solver_ref.executeStep(T_stop_final, early_stop=True)
    print(f"\u2713 Reference completed with {solver_ref.getNumSteps()} steps\n")

    filename = os.path.join(_THIS_DIR, "tenet_checkpoint.bin")

    # STEP 1: Run partial path and save checkpoint.
    solver1 = TENET_Solver(X, D, y, lambda2, normalize=True, intercept=True,
                           verbose=False)
    solver1.executeStep(T_stop_partial, early_stop=True)
    print(f"Checkpoint at partial stop: {solver1.getNumSteps()} steps")
    solver1.save(filename)
    print(f"\u2713 Checkpoint saved at '{filename}'")

    # STEP 2: Load from checkpoint and continue.
    # The Python binding hydrates an existing instance (cpp: static
    # load(filename, X, D)); the constructor re-binds the data views.
    solver2 = TENET_Solver(X, D, y, lambda2, normalize=True, intercept=True,
                           verbose=False)
    solver2.load(filename)
    print(f"Loaded from checkpoint: {solver2.getNumSteps()} steps")
    solver2.executeStep(T_stop_final, early_stop=True)

    print("\nCOMPARISON:")
    print(f"Reference steps: {solver_ref.getNumSteps()}")
    print(f"Reloaded steps:  {solver2.getNumSteps()}")
    print(f"RSS diff:   {abs(solver_ref.getRSS()[-1] - solver2.getRSS()[-1])}")
    print(f"R2 diff:    {abs(solver_ref.getR2()[-1] - solver2.getR2()[-1])}")
    paths_match = solver_ref.getActions() == solver2.getActions()
    print("\u2713 Paths match" if paths_match else "\u2717 Paths differ!")

    print("\nT-Elastic Net Diagnostics:")
    print(f"  Removals: {solver2.getNumRemovals()}")
    print(f"  Cycling ratio: {solver2.getCyclingRatio():.4f}")

    os.remove(filename)
    print("\u2713 Checkpoint file removed")
    print("\n")


# ==============================================================================
# Demo 4: Memory-Mapped Data
# ==============================================================================

def demo_memory_mapped(high_dim, rnd_coef, T_stop):
    print_section_header("Demo 4: T-Elastic Net with Memory-Mapped Data")

    n = 1000 if high_dim else 5000
    p = 5000 if high_dim else 1000  # cpp: 500,000 on disk; scaled for the port
    num_dummies = 10 * p
    snr = 1.0

    true_support = [4, 27, 149, 398, 420]
    true_coefs = [2.5, -0.4, -0.2, -0.8, 1.1] if rnd_coef else [1.0] * 5
    lambda2 = 0.5

    print_demo_config(n, p, num_dummies, T_stop, true_support, true_coefs, snr)

    print("Generating memory-mapped data...")
    X_file = os.path.join(_THIS_DIR, "demo_tenet_X.bin")
    D_file = os.path.join(_THIS_DIR, "demo_tenet_D.bin")
    y_file = os.path.join(_THIS_DIR, "demo_tenet_y.bin")

    X, D, y = gen_synthetic_data(n, p, num_dummies, true_support, true_coefs,
                                 snr, seed=42)

    try:
        # Write X/D/y to disk (cpp: SyntheticDataMapped) and drop the
        # in-memory copies; the solver reads the zero-copy mmap views.
        X_mmap = numpy_to_memmap(X_file, X)
        D_mmap = numpy_to_memmap(D_file, D)
        y_mmap = numpy_to_memmap(y_file, np.asfortranarray(y.reshape(-1, 1)))
        del X, D
        print("\u2713 Data generated on disk\n")

        print("Running T-Elastic Net on memory-mapped data...")
        solver = TENET_Solver(X_mmap.to_numpy(), D_mmap.to_numpy(), y_mmap.to_numpy(),
                              lambda2, normalize=True, intercept=True, verbose=True)
        solver.executeStep(T_stop, early_stop=True)
        print("\u2713 T-Elastic Net completed\n")

        print_selection(solver, true_support)
        print_selection_quality(solver, true_support)
        print("\nT-Elastic Net Diagnostics:")
        print(f"  Removals: {solver.getNumRemovals()}")
        print(f"  Cycling ratio: {solver.getCyclingRatio():.4f}")
    finally:
        print("\nCleaning up files...")
        for f in (X_file, D_file, y_file):
            if os.path.exists(f):
                os.remove(f)
        print("\u2713 All files removed\n")


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    print_section_header("T-Elastic Net Demo Suite")

    set_num_threads(6)
    print(f"Running with {get_max_threads()} threads\n")

    # Demo 1: Basic usage with internal normalization
    demo_early_stopping(high_dim=True, rnd_coef=False, T_stop=10)
    demo_early_stopping(high_dim=False, rnd_coef=False, T_stop=10)

    # Demo 2: External normalization
    demo_with_external_normalizer(high_dim=False, rnd_coef=False, T_stop=10)
    demo_with_external_normalizer(high_dim=True, rnd_coef=False, T_stop=10)

    # Demo 3: Serialization
    demo_serialization()

    # Demo 4: Memory-mapped
    demo_memory_mapped(high_dim=True, rnd_coef=False, T_stop=10)

    print_section_header("\u2713 All demos completed successfully")
