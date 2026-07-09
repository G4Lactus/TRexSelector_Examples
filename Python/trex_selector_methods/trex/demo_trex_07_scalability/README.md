# Demo 07: Scalability Benchmark (Python)

## Purpose

Benchmark the classical T-Rex selector's runtime and peak resident memory,
comparing **in-memory** execution against **chunked, memory-mapped
out-of-core** execution over an exponentially increasing `(n, p)` grid. Unlike
the C++ counterpart (`demo_trex_07_mc_sim_scalability`, an unbuilt placeholder),
this Python port is a working benchmark and mirrors the R Demo 07.

## DGP / Data

- Scenario grid targeting raw-`X` footprints (8 bytes/double) from ~1 GB to
  ~256 GB: `(5000, 25000)`, `(15000, 75000)`, `(30000, 150000)`,
  `(50000, 250000)`, `(80000, 400000)`
- `s = 10` (contiguous support `{0, ..., 9}`), SNR = 1.0, seed = 42, tFDR = 0.1
- Out-of-core `X` is generated column-chunk by column-chunk (`CHUNK_COLS = 1000`)
  straight into a `MemoryMappedMatrix`, so the full matrix never sits in RAM

## What It Computes

`run_scalability_benchmark()` runs each scenario twice:

- **In-Memory**: `_INMEM_CTRL` (TLARS, `use_memory_mapping = False`) on a
  fully materialized `dgp_gauss_snr` design.
- **Memory-Mapped**: `_MMAP_CTRL` (TLARS, `use_memory_mapping = True`,
  `lloop_strategy = "HCONCAT"`) on the chunked mmap-backed `X`.

Each run is executed in an **isolated subprocess** (spawn context) so its
`ru_maxrss` is its own peak and an out-of-memory scenario is reported as
`OOM/ERROR` and skipped rather than killing the whole demo. Execution time and
peak RSS are recorded per run.

## Imports

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `dgp_gauss_snr` and `trex_sim_common` modules. The package
imports as `trex_selector_neo` (`MemoryMappedMatrix`, `AccessMode`,
`_make_trex_ctrl_from_dict`).

## Output

Writes `d04_scalability_benchmark.csv` (columns:
`scenario, n, p, type, time_sec, mem_mb, status`) to this demo's own
`simulation_results/` folder, and prints the summary table.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex/demo_trex_07_scalability/demo_trex_07_scalability.py
```

Each benchmarked run executes in its own subprocess via Python `multiprocessing`
(spawn start method). Warning: the upper scenarios target very large footprints
(up to ~256 GB) and are expected to report OOM on typical hardware.

**Last updated**: 2026-07-08
