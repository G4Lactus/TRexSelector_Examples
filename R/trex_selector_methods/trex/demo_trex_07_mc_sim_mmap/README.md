# Demo 07: Monte Carlo Simulation with Memory-Mapped Data

## Purpose

Validate T-Rex performance under memory-mapped (mmap) I/O over many Monte Carlo trials,
extending the single-run patterns of Demo 06. Confirms that the storage medium is
transparent to the algorithm (FDR/TPR unchanged) and that the per-trial mmap file
lifecycle is exception-safe.

## Data-generating process

`dgp_gauss_snr` — i.i.d. Gaussian design, `y = X beta + eps` calibrated to SNR.
High-dimensional config: `n = 300`, `p = 1000`, `s = 10`, `tFDR = 0.1`. The true
support is drawn once and shared across all trials via
`make_support_random(s, p, seed = 24)`; coefficients fixed `beta_j = 1`. SNR sweep
`{0.1, 0.5, 1.0, 2.0, 5.0}`, TLARS only. `num_MC = 500` for both parts (the C++ demo
uses `num_MC = 10` for Part C and `500` for Part D). Control via
`make_mmap_trex_ctrl(solver = "TLARS", tloop_max_stagnant_steps = 7, opt_threshold = 0.75)`.

## Parts

- Part A (`run_part_a`, header "Part C") — MC with in-memory `X` + `use_memory_mapping = TRUE`:
  verifies the D-mmap + solver-serialization pipeline yields stable FDR/TPR over many runs.
- Part B (`run_part_b`, header "Part D") — MC with memory-mapped `X` + `use_memory_mapping = TRUE`:
  each trial writes its own temporary mmap-backed `X` via `execute_with_temp_mmap`; because
  `mclapply` forks, each trial runs in its own address space with no cross-trial file
  conflicts (mirrors the C++ RAII `MmapFileGuard`).

Both report averaged FDR, TPR, Avg L, Avg T per SNR (`track_L_T = TRUE`). MC trials run
in parallel over a hardcoded `num_cores <- 6L` (set at the top of the file).

## Running

```bash
Rscript R/trex_selector_methods/trex/demo_trex_07_mc_sim_mmap/demo_trex_07_mc_sim_mmap.R
```

Each part writes a `.txt` table and a tidy `.csv` (header `solver,metric,snr,value`) to
this demo's own `simulation_results/` folder, stems `d03_trex_mmap_demo_c_n300_p1000_sw7`
and `d03_trex_mmap_demo_d_n300_p1000_sw7`. The script sources `../../support_generators.R`,
`../../simulation_utils.R`, and `../trex_sim_utils.R`.

**Last updated**: 2026-07-08
