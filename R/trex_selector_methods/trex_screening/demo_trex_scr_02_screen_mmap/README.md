# Demo 02: Screen-TRex (Memory-Mapped)

## Purpose

Repeat Demo 01's Ordinary vs. Bootstrap-CI screening with **memory-mapped** storage, to validate that Screen-TRex's disk-backed workflow reproduces the in-memory behavior while enabling larger, RAM-constrained problems. Two mmap levels are exercised: `trex_control(use_memory_mapping = TRUE)` activates the internal dummy-matrix mmap + solver-serialization pipeline (D-mmap), and `convert_to_memory_mapped()` additionally moves the design matrix `X` itself to a memory-mapped file so it never has to reside fully in RAM. Mirrors [demo_trex_scr_02_screen_trex_mmap.cpp](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_02_screen_trex_mmap/).

## Data-generating process

`dgp_iid_snr` — i.i.d. Gaussian design with SNR-calibrated noise.

- Single runs (Parts A–B): `n = 300`, `p = 1000`, fixed 1-based support `{28, 150, 399, 5, 43}` (C++ 0-based `{27, 149, 398, 4, 42}`), `beta_j = 1.0`, `SNR = 5.0`, data seed 58, selector seed 42.
- MC sweep (Part C): `n = 300`, `p = 1000`, random `s = 10` support per trial, `SNR in {0.1, 0.5, 1.0, 2.0, 5.0}`, `num_MC = 40` (the C++ demo uses 500; reduced for a ~1-minute runtime). Results are statistically comparable to the C++ reference, not bit-identical — the two use different RNG streams.

Shared controls via `make_scr_trex_ctrl(use_memory_mapping = TRUE)`: `solver = "TLARS"`, `K = 20`.

## Methods compared

- **Screen-TRex Ordinary (mmap)** — `trex_screen_control(trex_method = "TREX")`
- **Screen-TRex Bootstrap-CI (mmap)** — same, with `use_bootstrap_CI = TRUE`

## Parts

- **Part A** — In-memory `X` + `use_memory_mapping = TRUE` (D-mmap only): one Ordinary and one Bootstrap-CI single run.
- **Part B** — Fully memory-mapped `X` (written to a temp mmap file via `convert_to_memory_mapped()`, unlinked on exit) + D-mmap: one Ordinary single run.
- **Part C** — Monte Carlo SNR sweep (Ordinary vs. Bootstrap-CI) with mmap enabled, reporting FDR / TPR / Estimated FDR.

## Running

```bash
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_02_screen_mmap/demo_trex_scr_02_screen_mmap.R
```

Runs from any working directory; the `run_part_*` flags at the top toggle individual parts; the script sources the shared `../trex_scr_sim_utils.R`. Part C writes `d02_screen_trex_mmap_mc_n300_p1000.txt` (summary table) and `.csv` (tidy long format) into the local [simulation_results/](simulation_results/) folder; Parts A–B print to the console only.

## Interpretation

- Best read side-by-side with Demo 01: DGP, controls, and methods are identical, so the FDR/TPR/estimated-FDR curves should match Demo 01 closely (up to Monte Carlo noise) — the only intended difference is where the dummy matrices (and, in Part B, `X`) are stored, not the statistical behavior of the selector. A material divergence at the same SNR points would suggest an mmap-specific issue rather than a genuine statistical effect.
- The recorded run (40 MC) shows the expected pattern: Ordinary FDR falls from 0.40 at SNR 0.1 to about 0.04 at SNR 5 with TPR reaching 1.00, while Bootstrap-CI is more conservative throughout (lower FDR, lower TPR — 0.85 at SNR 5).
- Parts A and B use the same data seed (58) and selector seed (42), so their Ordinary runs should agree — mmap-backing `X` is meant to be transparent to the selector.

**Last updated**: 2026-07-08
