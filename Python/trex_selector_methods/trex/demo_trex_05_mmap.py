# ==============================================================================
# demo_trex_05_mmap.py
# ==============================================================================
#
# Classical T-Rex Selector demo — memory-mapped workflows.
# Mirrors R/trex_selector_methods/trex/demo_trex_05_mmap.R and
# cpp/trex_selector_methods/trex/demo_trex_05_mmap.cpp.
#
# Demo content:
#
#  Part A: In-memory X + use_memory_mapping=True.
#          Setting use_memory_mapping=True activates both:
#            (1) Memory-mapped storage for internal dummy matrices D.
#            (2) Solver LARS-path serialization to disk between T-loop iterations.
#          X itself stays fully in RAM.
#          Exact mirror of C++ Demo A.
#
#  Part B: Fully memory-mapped pipeline.
#          X is written to a temporary binary file via numpy_to_memmap().
#          X_mmap.to_numpy() (a zero-copy OS-paged view) is passed to
#          TRexSelector; use_memory_mapping=True adds D-mmap + solver
#          serialization on top.
#          Mirrors C++ Demo B (SyntheticDataMapped) and R Part B
#          (execute_with_temp_mmap).
#
# Seed note:
#   Both parts use seed=58 for DGP and TRexSelector (single-run reproducibility,
#   mirrors R Part A/B).  This is intentional for a single-run demo; MC demos
#   use seed=-1 (hardware entropy) for valid FDR estimates.
#
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

import tempfile
import numpy as np
import trex_selector as ts
from trex_selector.utils import numpy_to_memmap

from dgp_gauss_snr import dgp_gauss_snr
from trex_sim_common import _make_trex_ctrl_from_dict
from trex_helpers import print_single_run_result, save_selection_csv

# ==============================================================================
# Shared demo configuration
# ==============================================================================
# Mirrors make_mmap_demo_config() in C++ and DEMO_CFG in R.
#
# Fixed true support (0-based) — note: different from demo_trex_01.
# C++ (0-based): {27, 149, 398, 420, 4}
# R   (1-based): {28, 150, 399, 421, 5}
_TRUE_SUPPORT = np.array([27, 149, 398, 420, 4], dtype=np.intp)
_TRUE_COEFS   = np.ones(5, dtype=np.float64)     # rnd_coef = False
_SNR          = 1.0
_TFDR         = 0.1
_SEED         = 58

_OUT_DIR = os.path.join(_THIS_DIR, "simulation_results")

# Shared control: mirrors make_mmap_trex_ctrl() in C++ and R
_MMAP_CTRL_DICT = dict(
    solver_name               = "TLARS",
    K                         = 20,
    max_dummy_multiplier      = 10,
    use_max_T_stop            = True,
    lloop_strategy            = "HCONCAT",
    tloop_stagnation_stop     = True,
    use_memory_mapping        = True,   # activates D-mmap + solver serialization
    lambda2                   = 0.1,
    rho_afs                   = 0.3,
    ncgmp_variant             = 0,
    parallel_rnd_experiments  = False,  # mirrors R use_openmp=FALSE
)


# ==============================================================================
# Shared print helper (mirrors R .print_and_save_d02)
# ==============================================================================

def _print_and_save(stem, dat, sel, tFDR):
    from trex_sim_common import compute_fdp, compute_tpp
    true_support = dat["true_support"]
    selected     = list(sel.selected_indices)
    n_sel        = len(selected)
    true_set     = set(true_support.tolist())
    n_tp         = len(set(selected) & true_set)
    n_fp         = n_sel - n_tp
    fdp_val      = compute_fdp(selected, true_support.tolist())
    tpp_val      = compute_tpp(selected, true_support.tolist())

    print("-" * 70)
    print(f"  Selected indices (0-based): {{{', '.join(str(x) for x in sorted(selected))}}}")
    print(f"  FDP = {fdp_val:.4f}  |  TPP = {tpp_val:.4f}  (target tFDR <= {tFDR:.2f})")
    print(f"  T_stop = {sel.T_stop},  L = {sel.L}")
    print("-" * 70 + "\n")

    save_selection_csv(
        _OUT_DIR, stem,
        p=dat["p"],
        phi_prime=sel.phi_prime,
        selected_0based=selected,
        true_support_0based=true_support.tolist(),
    )


# ==============================================================================
# Part A: In-memory X + use_memory_mapping=True
# ==============================================================================

def part_a_d_mmap_solver_serial():
    """
    Mirrors C++ demo_TRexSelector_d_mmap_solver_serial() and R trx_02_part_a().
    """
    def run_single(high_dim):
        n = 1000 if high_dim else 5000
        p = 5000 if high_dim else 1000
        dim_label = "High-dimensional (p > n)" if high_dim else "Low-dimensional (n > p)"

        print(f"\n{'=' * 70}")
        print(f"Part A: In-Memory X  +  use_memory_mapping=True")
        print(f"  Activates D-mmap and solver serialization inside TRexSelector.")
        print(f"  {dim_label}  (n = {n}, p = {p})\n")

        print("Generating synthetic data (in-memory) ...")
        dat = dgp_gauss_snr(n, p, _TRUE_SUPPORT, coefs=_TRUE_COEFS,
                            snr=_SNR, seed=_SEED)

        ctrl = _make_trex_ctrl_from_dict(_MMAP_CTRL_DICT)

        print("Running T-Rex Selector ...\n")
        sel = ts.TRexSelector(
            np.asfortranarray(dat["X"]),
            dat["y"],
            tFDR=_TFDR,
            trex_control=ctrl,
            seed=_SEED,
            verbose=True,
        )
        sel.select()

        stem = f"d02_trex_mmap_demo_a_n{n}_p{p}"
        _print_and_save(stem, dat, sel, _TFDR)

    run_single(high_dim=False)
    run_single(high_dim=True)


# ==============================================================================
# Part B: Memory-mapped X + use_memory_mapping=True
# ==============================================================================

def part_b_full_mmap():
    """
    Mirrors C++ demo_TRexSelector_full_mmap() and R trx_02_part_b().

    X is written to a temporary binary file; X_mmap.to_numpy() (a zero-copy
    Fortran-contiguous float64 view backed by the mmap buffer) is passed to
    TRexSelector.  The temp file is removed in a finally block — exception-safe,
    mirrors C++ RAII MmapFileGuard and R execute_with_temp_mmap / on.exit().
    """
    def run_single(high_dim):
        n = 1000 if high_dim else 5000
        p = 5000 if high_dim else 1000
        dim_label = "High-dimensional (p > n)" if high_dim else "Low-dimensional (n > p)"

        print(f"\n{'=' * 70}")
        print(f"Part B: Memory-Mapped X  +  use_memory_mapping=True")
        print(f"  Full mmap pipeline: X mmap + D mmap + solver serialization.")
        print(f"  {dim_label}  (n = {n}, p = {p})\n")

        print("Generating synthetic data ...")
        dat = dgp_gauss_snr(n, p, _TRUE_SUPPORT, coefs=_TRUE_COEFS,
                            snr=_SNR, seed=_SEED)

        ctrl = _make_trex_ctrl_from_dict(_MMAP_CTRL_DICT)

        # Create temporary backing file.  os.close(fd) releases the OS file
        # descriptor; numpy_to_memmap re-opens the path for writing.
        fd, mmap_path = tempfile.mkstemp(suffix=".bin")
        os.close(fd)

        try:
            print(f"Converting X to memory-mapped file ({mmap_path}) ...")
            X_mmap = numpy_to_memmap(mmap_path, dat["X"])

            print("Running T-Rex Selector ...\n")
            # X_mmap.to_numpy() → Fortran-contiguous float64 mmap view.
            # TRexSelector.__init__ calls np.asfortranarray(X) which returns
            # the view unchanged (already F-contiguous float64, zero copy).
            sel = ts.TRexSelector(
                X_mmap.to_numpy(),
                dat["y"],
                tFDR=_TFDR,
                trex_control=ctrl,
                seed=_SEED,
                verbose=True,
            )
            sel.select()
        finally:
            # Exception-safe cleanup.  On CPython the MemoryMappedMatrix is
            # reference-counted; it goes out of scope (and releases the mmap
            # mapping) before this finally block runs.
            try:
                os.unlink(mmap_path)
                print(f"[Info] Temp mmap file removed: {mmap_path}")
            except OSError:
                pass

        stem = f"d02_trex_mmap_demo_b_n{n}_p{p}"
        _print_and_save(stem, dat, sel, _TFDR)

    run_single(high_dim=False)
    run_single(high_dim=True)


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    part_a_d_mmap_solver_serial()
    part_b_full_mmap()
    print("\ndemo_trex_05_mmap.py complete.")
