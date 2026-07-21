# ==============================================================================
# demo_trex_08_mc_sim_scalability.py
# ==============================================================================
#
# Scalability benchmark for the T-Rex Selector: runtime and memory usage
# evaluated as n and p scale, across three base solvers (TLARS, TOMP, TAFS).
# Mirrors cpp/trex_selector_methods/trex/demo_trex_08_mc_sim_scalability.cpp.
#
# Runtime and memory benchmark on the fully memory-mapped pipeline: every MC
# trial generates X directly into an mmap backing file (chunked column
# generation via MemoryMappedMatrix, so X is never materialized in RAM), and
# use_memory_mapping=True additionally memory-maps the dummy matrices D and
# serializes solver warm-start state between T-loop iterations.
#
#  Scenarios (SCENARIOS_BENCHMARK):
#    A — n =   300, p =   1,000  (baseline)
#    B — n = 1,000, p =  10,000  (moderate scale)
#    C — n = 5,000, p = 100,000  (large sample x high-dimensional)
#
#  2D sweep: every scenario runs over the SNR grid {0.5, 1.0, 2.0} for the
#  base solvers TLARS, TOMP, and TAFS (rho_afs = 0.3), with the demo-03 DGP
#  (variable support of cardinality s = 10 redrawn per trial, unit
#  coefficients, Normal predictors/noise).
#
#  Per grid point (scenario x solver x SNR; num_MC SEQUENTIAL trials so
#  wall-clock and RSS readings are undistorted — this demo intentionally does
#  NOT use the ProcessPoolExecutor MC runners): timing breakdown (DGP,
#  select(), per-T-iteration, total), memory (post-select RSS, physical
#  footprint, process peak RSS, X/y backing sizes), and FDR / TPR / Avg L /
#  Avg T.
#
#  Results go to simulation_results/data/ as an aligned table (.txt) and a
#  tidy CSV, re-saved after each completed scenario — a crash on a larger
#  scenario preserves the finished ones.
#
# Python adaptations (noted, cpp is the reference):
#  - y is kept in RAM (n doubles — negligible); the "X+y MiB" metric counts
#    the X backing file plus y's in-memory bytes.
#  - X is filled column-chunk by column-chunk into the mmap view, so the RNG
#    stream mapping (and hence the exact data) differs from dgp_gauss_snr and
#    from the C++ SyntheticDataMapped — irrelevant for a benchmark.
#  - Memory hooks: RSS via `ps` (macOS) / /proc/self/statm (Linux); physical
#    footprint via task_info(TASK_VM_INFO).phys_footprint through ctypes
#    (macOS) / VmRSS + VmSwap (Linux); peak RSS via resource.getrusage.
#
# NOTE: Scenario C is heavy (minutes per trial, tens of GiB of
# swap/compression) — mirror of the C++ demo's runtime warning.
#
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import ctypes
import resource
import subprocess
import tempfile
import time

import numpy as np
import pandas as pd

import trex_selector_neo as tsn
from trex_selector_neo.utils import AccessMode, MemoryMappedMatrix

from trex_helpers import make_dual_out
from trex_sim_common import _make_trex_ctrl_from_dict, compute_fdp, compute_tpp
from support_generators import make_support_random

# ==============================================================================
# Configuration
# ==============================================================================

_OUT_DIR = os.path.join(_THIS_DIR, "simulation_results", "data")

_S_TRUE     = 10     # cardinality of the true support (redrawn per trial)
_TFDR       = 0.1
_SNR_VALUES = [0.5, 1.0, 2.0]
_CHUNK_COLS = 1000   # column-chunk size for the out-of-core X generation

# Scenario tuples: (label, n, p, num_MC)  — mirrors make_benchmark_scenarios()
SCENARIOS_BENCHMARK = [
    ("A", 300,  1_000,   10),
    ("B", 1_000, 10_000, 10),
    ("C", 5_000, 100_000, 10),
]

# Tiny end-to-end smoke-test scenario (pipeline verification only)
SCENARIOS_SMOKE = [
    ("S", 150, 500, 2),
]

# Base solvers to compare (demo 05 configuration; rho_afs only used by TAFS)
_SOLVERS = [
    dict(solver_name="TLARS", name="TLARS", lambda2=0.1, rho_afs=0.3, ncgmp_variant=0),
    dict(solver_name="TOMP",  name="TOMP",  lambda2=0.1, rho_afs=0.3, ncgmp_variant=0),
    dict(solver_name="TAFS",  name="TAFS",  lambda2=0.1, rho_afs=0.3, ncgmp_variant=0),
]


def _make_scalability_ctrl_dict(solver_desc):
    """T-Rex control (base config of demo 02/03) for one grid point.

    Mirrors make_scalability_trex_ctrl(): always-on use_memory_mapping
    (D mmap + solver serialization) on top of the mmap-generated X.
    parallel_rnd_experiments keeps its default (True) — trials are
    sequential, but the random experiments inside select() may run in
    parallel, mirroring the C++ OpenMP setup.
    """
    return dict(
        solver_name              = solver_desc["solver_name"],
        lambda2                  = solver_desc["lambda2"],
        rho_afs                  = solver_desc["rho_afs"],
        ncgmp_variant            = solver_desc["ncgmp_variant"],
        K                        = 20,
        max_dummy_multiplier     = 10,
        use_max_T_stop           = True,
        lloop_strategy           = "STANDARD",
        tloop_stagnation_stop    = True,
        tloop_max_stagnant_steps = 5,
        use_memory_mapping       = True,
    )


# ==============================================================================
# Memory profiling hooks (mirror current_rss_mib / current_footprint_mib /
# peak_rss_mib in the C++ demo)
# ==============================================================================

def current_rss_mib():
    """Current resident set size (RSS) of this process in MiB (0.0 if unavailable)."""
    if sys.platform == "darwin":
        try:
            out = subprocess.run(["ps", "-o", "rss=", "-p", str(os.getpid())],
                                 capture_output=True, text=True, check=True)
            return int(out.stdout.strip()) / 1024.0        # ps reports KiB
        except Exception:
            return 0.0
    try:
        with open("/proc/self/statm") as fh:
            resident_pages = int(fh.read().split()[1])
        return resident_pages * os.sysconf("SC_PAGESIZE") / (1024.0 * 1024.0)
    except Exception:
        return 0.0


class _TaskVMInfo(ctypes.Structure):
    """mach task_vm_info through the phys_footprint field (TASK_VM_INFO rev1)."""
    _fields_ = [
        ("virtual_size",                ctypes.c_uint64),
        ("region_count",                ctypes.c_int32),
        ("page_size",                   ctypes.c_int32),
        ("resident_size",               ctypes.c_uint64),
        ("resident_size_peak",          ctypes.c_uint64),
        ("device",                      ctypes.c_uint64),
        ("device_peak",                 ctypes.c_uint64),
        ("internal",                    ctypes.c_uint64),
        ("internal_peak",               ctypes.c_uint64),
        ("external",                    ctypes.c_uint64),
        ("external_peak",               ctypes.c_uint64),
        ("reusable",                    ctypes.c_uint64),
        ("reusable_peak",               ctypes.c_uint64),
        ("purgeable_volatile_pmap",     ctypes.c_uint64),
        ("purgeable_volatile_resident", ctypes.c_uint64),
        ("purgeable_volatile_virtual",  ctypes.c_uint64),
        ("compressed",                  ctypes.c_uint64),
        ("compressed_peak",             ctypes.c_uint64),
        ("compressed_lifetime",         ctypes.c_uint64),
        ("phys_footprint",              ctypes.c_uint64),
    ]


def current_footprint_mib():
    """
    Physical footprint of this process in MiB: anonymous heap including
    compressed / swapped-out pages (Activity Monitor's "Memory" on macOS).
    Linux approximation: VmRSS + VmSwap. Falls back to RSS on failure.
    """
    if sys.platform == "darwin":
        try:
            libc = ctypes.CDLL(None)
            info  = _TaskVMInfo()
            count = ctypes.c_uint32(ctypes.sizeof(_TaskVMInfo) // 4)
            TASK_VM_INFO = 22
            kr = libc.task_info(libc.mach_task_self(), TASK_VM_INFO,
                                ctypes.byref(info), ctypes.byref(count))
            if kr == 0 and info.phys_footprint > 0:
                return info.phys_footprint / (1024.0 * 1024.0)
        except Exception:
            pass
        return current_rss_mib()
    try:
        kib = 0.0
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith(("VmRSS:", "VmSwap:")):
                    kib += float(line.split()[1])
        return kib / 1024.0
    except Exception:
        return current_rss_mib()


def peak_rss_mib():
    """Peak RSS of this process in MiB (monotone over the process lifetime)."""
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss is bytes on macOS, KiB on Linux
    return r / (1024.0 * 1024.0) if sys.platform == "darwin" else r / 1024.0


# ==============================================================================
# Memory-mapped DGP (mirrors SyntheticDataMapped: X generated on disk)
# ==============================================================================

def dgp_gauss_snr_mmap(n, p, support, coefs, snr, seed, x_path,
                       chunk_cols=_CHUNK_COLS):
    """
    Generate X column-chunk by column-chunk into an mmap backing file and an
    SNR-calibrated y (kept in RAM; n doubles are negligible). X never sits
    fully in RAM; working memory is O(n * chunk_cols).

    Returns (X_mmap, y) with X_mmap a MemoryMappedMatrix over x_path.
    """
    rng  = np.random.default_rng(seed)
    beta = np.zeros(p, dtype=np.float64)
    beta[np.asarray(support, dtype=np.intp)] = np.asarray(coefs, dtype=np.float64)

    X = MemoryMappedMatrix(x_path, n, p, AccessMode.ReadWrite)
    V = X.to_numpy()                     # zero-copy F-contiguous view
    signal = np.zeros(n)
    for s0 in range(0, p, chunk_cols):
        e0    = min(s0 + chunk_cols, p)
        chunk = rng.standard_normal((n, e0 - s0))
        V[:, s0:e0] = chunk
        signal += chunk @ beta[s0:e0]

    noise_sigma = np.sqrt(np.var(signal, ddof=1) / snr)
    y = signal + rng.standard_normal(n) * noise_sigma
    return X, y


# ==============================================================================
# Single-trial runner (sequential; timing + memory instrumented)
# ==============================================================================

def run_scalability_trial(n, p, ctrl_dict, tFDR, snr, seed, support, coefs,
                          mmap_filestem):
    """
    One MC trial on the fully memory-mapped pipeline. The X backing file is
    removed in a finally block (exception-safe; mirrors the C++ RAII
    MmapFileGuard). Returns a dict of per-trial measurements.
    """
    m = {}
    x_path = f"{mmap_filestem}_X.dat"

    try:
        t0 = time.perf_counter()
        X_mmap, y = dgp_gauss_snr_mmap(n, p, support, coefs, snr, seed, x_path)
        m["dgp_s"] = time.perf_counter() - t0

        m["data_mib"] = (os.path.getsize(x_path) + y.nbytes) / (1024.0 * 1024.0)

        ctrl = _make_trex_ctrl_from_dict(ctrl_dict)

        t1 = time.perf_counter()
        # X_mmap.to_numpy() → zero-copy F-contiguous float64 mmap view;
        # seed=-1 (hardware entropy) as in all MC demos.
        sel = tsn.TRexSelector(X_mmap.to_numpy(), y, tFDR=tFDR,
                               trex_control=ctrl, seed=-1, verbose=False)
        sel.select()
        m["select_s"] = time.perf_counter() - t1

        m["rss_mib"] = current_rss_mib()
        m["fp_mib"]  = current_footprint_mib()

        selected = list(sel.selected_indices)
        support_list = list(support)
        m["fdp"] = compute_fdp(selected, support_list)
        m["tpp"] = compute_tpp(selected, support_list)
        m["L"]   = float(sel.L)
        m["T"]   = float(sel.T_stop)

        # Drop the selector and mmap view before the backing file is removed
        # (CPython ref-counting makes this deterministic).
        del sel, X_mmap
    finally:
        try:
            os.unlink(x_path)
        except OSError:
            pass
    return m


# ==============================================================================
# Grid-point loop (scenario x solver x SNR): sequential MC trials
# ==============================================================================

def run_scalability_grid_point(scenario, solver_desc, tFDR, snr, base_seed,
                               rnd_coef, cardinality_true_support):
    label, n, p, num_MC = scenario
    ctrl_dict = _make_scalability_ctrl_dict(solver_desc)

    mmap_filestem = os.path.join(tempfile.gettempdir(),
                                 f"demo08_scalability_{label}")

    print(f"  [{label}/{solver_desc['name']}/SNR={snr:.1f}] "
          f"n={n}, p={p} — running {num_MC} sequential MC trials (mmap) ...",
          flush=True)

    agg = dict(avg_fdr=0.0, avg_tpr=0.0, avg_L=0.0, avg_T=0.0,
               avg_dgp_s=0.0, avg_select_s=0.0, avg_titer_s=0.0,
               avg_total_s=0.0, max_rss_mib=0.0, max_fp_mib=0.0,
               peak_rss_mib=0.0, data_mib=0.0)

    for mc in range(num_MC):
        trial_seed = base_seed + mc

        # Variable support (demo 03 pattern; make_support_random adds the
        # +500_000 offset internally, mirroring C++ draw_support)
        support = make_support_random(cardinality_true_support, p,
                                      seed=trial_seed)
        if rnd_coef:
            coef_rng = np.random.default_rng(int(trial_seed) + 600_000)
            coefs = coef_rng.standard_normal(cardinality_true_support)
        else:
            coefs = np.ones(cardinality_true_support)

        m = run_scalability_trial(n, p, ctrl_dict, tFDR, snr, trial_seed,
                                  support, coefs, mmap_filestem)

        agg["avg_fdr"]      += m["fdp"]
        agg["avg_tpr"]      += m["tpp"]
        agg["avg_L"]        += m["L"]
        agg["avg_T"]        += m["T"]
        agg["avg_dgp_s"]    += m["dgp_s"]
        agg["avg_select_s"] += m["select_s"]
        agg["avg_titer_s"]  += m["select_s"] / max(m["T"], 1.0)
        agg["avg_total_s"]  += m["dgp_s"] + m["select_s"]
        agg["max_rss_mib"]   = max(agg["max_rss_mib"], m["rss_mib"])
        agg["max_fp_mib"]    = max(agg["max_fp_mib"], m["fp_mib"])
        agg["data_mib"]      = m["data_mib"]

        print(f"    trial {mc + 1}/{num_MC}: dgp={m['dgp_s']:.2f}s  "
              f"select={m['select_s']:.2f}s  FDP={m['fdp']:.3f}  "
              f"TPP={m['tpp']:.3f}  L={m['L']:.0f}  T={m['T']:.0f}  "
              f"RSS={m['rss_mib']:.1f} MiB  FP={m['fp_mib']:.1f} MiB",
              flush=True)

    for key in ("avg_fdr", "avg_tpr", "avg_L", "avg_T",
                "avg_dgp_s", "avg_select_s", "avg_titer_s", "avg_total_s"):
        agg[key] /= float(num_MC)
    agg["peak_rss_mib"] = peak_rss_mib()

    print(f"  [{label}/{solver_desc['name']}/SNR={snr:.1f}] done. "
          f"avg select={agg['avg_select_s']:.2f}s  FDR={agg['avg_fdr']:.3f}  "
          f"TPR={agg['avg_tpr']:.3f}\n", flush=True)
    return agg


# ==============================================================================
# Output: aligned table (console + .txt) and tidy CSV
# ==============================================================================

# Metric rows reported per solver (table row order and CSV order):
# (label, table decimals, aggregate key)
_METRICS = [
    ("FDR",       4, "avg_fdr"),
    ("TPR",       4, "avg_tpr"),
    ("Avg L",     2, "avg_L"),
    ("Avg T",     2, "avg_T"),
    ("DGP s",     3, "avg_dgp_s"),
    ("Select s",  3, "avg_select_s"),
    ("s/T-iter",  4, "avg_titer_s"),
    ("Total s",   3, "avg_total_s"),
    ("RSS MiB",   1, "max_rss_mib"),
    ("FootprMiB", 1, "max_fp_mib"),
    ("PeakRSS",   1, "peak_rss_mib"),
    ("X+y MiB",   1, "data_mib"),
]


def save_and_print_scalability_results(file_stem, tFDR, snr_values, scenarios,
                                       solvers, results):
    """
    Print averaged scalability results as an aligned table (console + .txt)
    and write a tidy long-format CSV (demo 02--07 table style; grid columns
    are scenario x SNR combinations). Mirrors the C++ function of the same
    name; `results` maps solver name -> flat list of grid-point dicts indexed
    as sc_idx * len(snr_values) + snr_idx.
    """
    S = len(snr_values)
    os.makedirs(_OUT_DIR, exist_ok=True)
    txt_path = os.path.join(_OUT_DIR, f"{file_stem}.txt")

    with open(txt_path, "w") as fh:
        out = make_dual_out(fh)

        # Header + scenario info block
        snr_str = ", ".join(f"{x:.1f}" for x in snr_values)
        out("\n")
        out("=" * 70 + "\n")
        out(f"=== T-Rex Scalability Results (tFDR = {tFDR:.2f}, "
            f"SNR grid = {{{snr_str}}}) ===\n")
        out("=" * 70 + "\n\n")
        out("All scenarios: fully memory-mapped pipeline "
            "(X mmap + D mmap + solver serialization)\n\n")
        out(f"{'Scenario':<10}{'n':<9}{'p':<10}{'num_MC':<8}\n")
        out("-" * 37 + "\n")
        for label, n, p, num_MC in scenarios:
            out(f"{label:<10}{n:<9}{p:<10}{num_MC:<8}\n")
        out("\n")

        # Column header (one column per scenario x SNR, labelled "<sc>/<snr>")
        indent_w, metric_w, col_w = 2, 23, 14
        sep_w = indent_w + metric_w + col_w * len(scenarios) * S
        hdr = " " * indent_w + f"{'Scenario/SNR':<{metric_w}}"
        for label, _, _, _ in scenarios:
            for snr in snr_values:
                hdr += f"{f'{label}/{snr:.1f}':>{col_w}}"
        out(hdr + "\n" + "-" * sep_w + "\n")

        # Data rows: per solver, one row per metric
        for solver in solvers:
            grid = results[solver["name"]]
            out("\n" + solver["name"] + "\n")
            for name, prec, key in _METRICS:
                row = " " * indent_w + f"{name:<{metric_w}}"
                for i in range(len(scenarios) * S):
                    row += f"{grid[i][key]:>{col_w}.{prec}f}"
                out(row + "\n")
        out("\n")

    # Tidy long-format CSV (scenario, solver, n, p, num_mc, snr, metric, value)
    rows = []
    for solver in solvers:
        grid = results[solver["name"]]
        for i, (label, n, p, num_MC) in enumerate(scenarios):
            for si, snr in enumerate(snr_values):
                for name, _, key in _METRICS:
                    rows.append(dict(scenario=label, solver=solver["name"],
                                     n=n, p=p, num_mc=num_MC,
                                     snr=f"{snr:.6f}", metric=name,
                                     value=round(grid[i * S + si][key], 6)))
    csv_path = os.path.join(_OUT_DIR, f"{file_stem}.csv")
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    print(f"[Info] CSV results saved to:             {csv_path}")
    print(f"[Info] Results successfully saved to: {txt_path}\n")


# ==============================================================================
# Scalability benchmark driver
# ==============================================================================

def demo_trex_scalability(scenarios, file_stem, rnd_coef=False):
    print("\n" + "=" * 70)
    print("Demo: T-Rex Selector Scalability Benchmark  |  TLARS / TOMP / TAFS")
    print("=" * 70 + "\n")

    print("Fully memory-mapped pipeline in all scenarios "
          "(X mmap + D mmap + solver serialization)")
    print("Scenarios: " + "  ".join(
        f"{lbl} (n={n}, p={p}, num_MC={mc})" for lbl, n, p, mc in scenarios))
    print(f"SNR grid = {{{', '.join(f'{x:.1f}' for x in _SNR_VALUES)}}},  "
          f"s = {_S_TRUE},  tFDR = {_TFDR:.2f}\n")

    # Results: solver name -> grid points flattened as scenario x SNR
    S = len(_SNR_VALUES)
    empty = {key: 0.0 for _, _, key in _METRICS}
    results = {sv["name"]: [dict(empty) for _ in range(len(scenarios) * S)]
               for sv in _SOLVERS}

    # Scenario-outer loop: every solver finishes the current (smaller)
    # scenario before the next (larger) one starts; results are re-saved
    # after each completed scenario.
    for sc_idx, scenario in enumerate(scenarios):
        label, n, p, _ = scenario
        print("=" * 50)
        print(f"Scenario {label}: n={n}, p={p}")
        print("=" * 50)
        X_gib = n * p * 8.0 / (1024.0 ** 3)
        print(f"  Expected X backing file: {X_gib:.2f} GiB "
              f"(mmap D grows up to 10x that)\n", flush=True)

        # 2D sweep solver x SNR: the base seed is unique per (scenario, SNR)
        # and shared across solvers -> every solver sees identical data.
        for solver_desc in _SOLVERS:
            for snr_idx, snr in enumerate(_SNR_VALUES):
                base_seed = 24 + (sc_idx * S + snr_idx) * 1000
                results[solver_desc["name"]][sc_idx * S + snr_idx] = \
                    run_scalability_grid_point(scenario, solver_desc, _TFDR,
                                               snr, base_seed, rnd_coef,
                                               _S_TRUE)

        # Incremental save over completed scenarios
        save_and_print_scalability_results(file_stem, _TFDR, _SNR_VALUES,
                                           scenarios[:sc_idx + 1], _SOLVERS,
                                           results)


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    # Smoke test: tiny sizes, verifies the full pipeline end-to-end
    if False:
        demo_trex_scalability(SCENARIOS_SMOKE,
                              "demo_trex_08_scalability_smoke")

    # Scalability benchmark: scenarios A / B / C
    if True:
        demo_trex_scalability(SCENARIOS_BENCHMARK,
                              "demo_trex_08_scalability_results")

    print("\ndemo_trex_08_mc_sim_scalability.py complete.")
