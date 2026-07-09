# ==============================================================================
# demo_trex_07_scalability.py
# ==============================================================================
#
# Part 7 of the classical T-Rex demos: a scalability benchmark comparing
# in-memory execution against chunked, memory-mapped out-of-core execution
# over an exponentially increasing (n, p) grid.
#
# Mirrors R/trex_selector_methods/trex/demo_trex_07_scalability.R.
#
# Demo content:
#   - A grid of (n, p) scenarios targeting raw-X footprints from 1 GB to 256 GB.
#   - Each scenario runs the classical T-Rex once in-memory and once out-of-core
#     (X generated column-chunk by column-chunk straight into an mmap file, so
#     the full matrix never sits in RAM).
#   - Execution time and peak resident memory are recorded per run.
#   - Results are written to a CSV and printed.
#
# Peak-memory note: R resets its GC max-used counter per run; Python has no
# resettable analog that also captures numpy / C++ allocations. Each run is
# therefore executed in an isolated subprocess whose ru_maxrss IS its own peak.
# The subprocess also lets an out-of-memory scenario be reported (OOM) and
# skipped, rather than killing the whole demo.
#
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import csv
import time
import resource
import tempfile
import multiprocessing

import numpy as np
import trex_selector_neo as ts
from trex_selector_neo.utils import MemoryMappedMatrix, AccessMode

from dgp_gauss_snr import dgp_gauss_snr
from trex_sim_common import _make_trex_ctrl_from_dict

# ==============================================================================
# Configuration
# ==============================================================================
# Exponential (n, p) grid targeting raw-X footprints (8 bytes per double).
SCENARIOS = [
    (5000,  25000,  "1 GB"),
    (15000, 75000,  "9 GB"),
    (30000, 150000, "36 GB"),
    (50000, 250000, "100 GB"),
    (80000, 400000, "256 GB"),
]

S_TRUE     = 10
SNR        = 1.0
SEED       = 42
TFDR       = 0.1
CHUNK_COLS = 1000

_OUT_DIR = os.path.join(_THIS_DIR, "simulation_results")

# In-memory: plain TLARS. Out-of-core: mmap control (mirrors make_mmap_trex_ctrl).
_INMEM_CTRL = dict(
    solver_name              = "TLARS",
    use_memory_mapping       = False,
    parallel_rnd_experiments = False,
)
_MMAP_CTRL = dict(
    solver_name              = "TLARS",
    K                        = 20,
    max_dummy_multiplier     = 10,
    use_max_T_stop           = True,
    lloop_strategy           = "HCONCAT",
    tloop_stagnation_stop    = True,
    use_memory_mapping       = True,
    parallel_rnd_experiments = False,
)


# ==============================================================================
# Workers (run in an isolated subprocess)
# ==============================================================================

def _peak_rss_mb():
    # ru_maxrss is bytes on macOS, kilobytes on Linux.
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 * 1024) if sys.platform == "darwin" else r / 1024


def _chunked_gauss_snr(n, p, support, x_path, chunk_cols, snr, seed):
    """Generate X column-chunk by column-chunk into an mmap file; return the
    mmap matrix and an SNR-calibrated y. X is never fully materialised in RAM."""
    rng  = np.random.default_rng(seed)
    beta = np.zeros(p)
    beta[support] = 1.0

    X = MemoryMappedMatrix(x_path, n, p, AccessMode.ReadWrite)
    V = X.to_numpy()                 # zero-copy F-contiguous view
    signal = np.zeros(n)
    for s in range(0, p, chunk_cols):
        e     = min(s + chunk_cols, p)
        chunk = rng.standard_normal((n, e - s))
        V[:, s:e] = chunk
        signal += chunk @ beta[s:e]

    noise_sigma = np.sqrt(np.var(signal, ddof=1) / snr)
    y = signal + rng.standard_normal(n) * noise_sigma
    return X, y


def _bench_child(payload, q):
    mode, n, p, ctrl_dict = payload
    support = np.arange(S_TRUE, dtype=np.intp)
    ctrl    = _make_trex_ctrl_from_dict(ctrl_dict)

    t0 = time.perf_counter()
    if mode == "inmem":
        dat = dgp_gauss_snr(n, p, support, snr=SNR, seed=SEED)
        sel = ts.TRexSelector(np.asfortranarray(dat["X"]), dat["y"],
                              tFDR=TFDR, trex_control=ctrl, seed=SEED, verbose=False)
        sel.select()
    else:
        fd, x_path = tempfile.mkstemp(suffix=".dat")
        os.close(fd)
        os.unlink(x_path)            # let MemoryMappedMatrix create/size it
        try:
            X, y = _chunked_gauss_snr(n, p, support, x_path, CHUNK_COLS, SNR, SEED)
            sel  = ts.TRexSelector(X.to_numpy(), y, tFDR=TFDR,
                                   trex_control=ctrl, seed=SEED, verbose=False)
            sel.select()
        finally:
            if os.path.exists(x_path):
                os.unlink(x_path)

    q.put((time.perf_counter() - t0, _peak_rss_mb(), "OK"))


def _run_isolated(mode, n, p, ctrl_dict):
    ctx = multiprocessing.get_context("spawn")
    q   = ctx.Queue()
    pr  = ctx.Process(target=_bench_child, args=((mode, n, p, ctrl_dict), q))
    pr.start()
    pr.join()
    if pr.exitcode != 0 or q.empty():
        return (float("nan"), float("nan"), "OOM/ERROR")
    return q.get()


# ==============================================================================
# Benchmark loop
# ==============================================================================

def run_scalability_benchmark():
    print("\n" + "=" * 70)
    print("Scalability Benchmark: In-Memory vs. Memory-Mapped")
    print("=" * 70 + "\n")

    rows = []
    for i, (n, p, lbl) in enumerate(SCENARIOS, 1):
        print(f"Scenario {i}: n = {n}, p = {p} (~{lbl} of raw data)")
        print("-" * 50)

        print("  [In-Memory]")
        t, m, st = _run_isolated("inmem", n, p, _INMEM_CTRL)
        print(f"    Status: {st} | Time: {t:.1f} s | Max RAM: {m:.0f} MB")
        rows.append((lbl, n, p, "In-Memory", t, m, st))

        print("  [Memory-Mapped]")
        t, m, st = _run_isolated("mmap", n, p, _MMAP_CTRL)
        print(f"    Status: {st} | Time: {t:.1f} s | Max RAM: {m:.0f} MB\n")
        rows.append((lbl, n, p, "Memory-Mapped", t, m, st))

    os.makedirs(_OUT_DIR, exist_ok=True)
    csv_path = os.path.join(_OUT_DIR, "d04_scalability_benchmark.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["scenario", "n", "p", "type", "time_sec", "mem_mb", "status"])
        w.writerows(rows)
    print(f"[Info] Benchmark finished. Results saved to: {csv_path}\n")

    hdr = ("scenario", "n", "p", "type", "time_sec", "mem_mb", "status")
    print(f"{hdr[0]:<10}{hdr[3]:<16}{hdr[4]:>10}{hdr[5]:>12}  {hdr[6]}")
    for r in rows:
        print(f"{r[0]:<10}{r[3]:<16}{r[4]:>10.1f}{r[5]:>12.0f}  {r[6]}")


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    run_scalability_benchmark()
    print("\ndemo_trex_07_scalability.py complete.")
