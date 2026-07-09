# Demo 07: Scalability Benchmark — In-Memory vs. Memory-Mapped

## Purpose

Benchmark the scalability of the T-Rex selector by comparing **in-memory** execution
against **chunked, memory-mapped out-of-core** execution over an exponentially growing
`(n, p)` grid, measuring wall-clock time and peak memory for each.

> Note: unlike the other demos, the C++ counterpart
> (`cpp/trex_selector_methods/trex/demo_trex_07_mc_sim_scalability/`) is an unimplemented
> placeholder (empty source, excluded from the build). This R port is the fully
> implemented version.

## Data-generating process

TLARS solver, `s_true = 10` (support `seq_len(s)`), `SNR = 1.0`, `seed = 42`,
`tFDR = 0.1`. In-memory data via `dgp_gauss_snr`; out-of-core data via
`dgp_chunked_gauss_snr(..., chunk_cols = 1000)` written to a temporary `.dat` file
(removed on exit). Exponential `(n, p)` grid `SCENARIOS`, sized to target these raw-`X`
footprints (8 bytes/double):

| Scenario | n | p | approx. X size |
|----------|------|--------|------|
| 1 | 5000 | 25000 | 1 GB |
| 2 | 15000 | 75000 | 9 GB |
| 3 | 30000 | 150000 | 36 GB |
| 4 | 50000 | 250000 | 100 GB |
| 5 | 80000 | 400000 | 256 GB |

## What it computes

- `run_scalability_benchmark` (`run_part_7`) iterates the scenarios; per scenario it runs
  `run_in_memory` (control `trex_control(solver = "TLARS")`) and `run_mmap` (control
  `make_mmap_trex_ctrl(solver = "TLARS")`), each wrapped in `tryCatch` so an OOM/error is
  recorded as a failed status rather than aborting the sweep.
- Records `time_sec` and peak `mem_mb` (max used from `gc()`), plus a `status`, per
  approach.

This is a resource benchmark, not an FDR/TPR study; the larger scenarios require
substantial disk/RAM and may be skipped as failures on modest hardware.

## Running

```bash
Rscript R/trex_selector_methods/trex/demo_trex_07_scalability/demo_trex_07_scalability.R
```

Results are written to this demo's own `simulation_results/` folder as
`d04_scalability_benchmark.csv` (columns `scenario,n,p,type,time_sec,mem_mb,status`) and
printed to the console. The script sources `../trex_sim_utils.R`.

**Last updated**: 2026-07-08
