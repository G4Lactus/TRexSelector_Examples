# Demo 01: Screen-TRex (In-Memory)

## Purpose

Basic in-memory Screen-TRex variable screening: two single-run walkthroughs of the **Ordinary** (`Phi_j > 0.5`) and **Bootstrap-CI** thresholding rules on the same fixed-support dataset, followed by a Monte Carlo SNR sweep comparing the two on FDR, TPR, and the procedure's self-reported Estimated FDR. In-memory variant — see Demo 02 for the memory-mapped counterpart. Mirrors `cpp/trex_selector_methods/trex_screening/demo_trex_scr_01_screen_trex/`.

## Data-generating process

`dgp_iid_snr` — i.i.d. standard-Gaussian design with SNR-calibrated noise, `sigma = sqrt(||X beta||^2 / (n * SNR))`.

- **Single runs (Parts 1–2)**: `n = 300`, `p = 1000` (high-dimensional, `p > n`), fixed 1-based true support `{28, 150, 399, 5, 43}` (`s = 5`; the C++ 0-based `{27, 149, 398, 4, 42}`), coefficients `beta_j = 1`, `SNR = 5.0`, data seed 58, selector seed 42.
- **MC sweep (Part 3)**: `n = 300`, `p = 1000`, random support of size `s = 10` drawn fresh per trial, `SNR in {0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`, `num_MC = 60` (the C++ demo uses 500; reduced for a ~1-minute runtime).

## Methods compared

- **Screen-TRex Ordinary** — `trex_screen_control(trex_method = "TREX", use_bootstrap_CI = FALSE)`
- **Screen-TRex Bootstrap-CI** — same, with `use_bootstrap_CI = TRUE` (library defaults for `R_boot` / `ci_grid_step`)

Shared solver controls via `make_scr_trex_ctrl()`: `solver = "TLARS"`, `K = 20`, `max_dummy_multiplier = 10`, `use_max_T_stop = TRUE`, `lloop_strategy = "STANDARD"`, `use_memory_mapping = FALSE`.

## Parts

- **Part 1** — Single-run Ordinary Screen-TRex: selected indices (1-based), FDP/TPP against the known support, Estimated FDR.
- **Part 2** — Single-run Bootstrap-CI Screen-TRex on the same data.
- **Part 3** — Monte Carlo SNR sweep, Ordinary vs. Bootstrap-CI, averaging FDR / TPR / Estimated FDR over 60 runs per (method, SNR) point.

## Running

```bash
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_01_screen/demo_trex_scr_01_screen.R
```

Runs from any working directory; the `run_part_*` flags at the top toggle individual parts; the script sources the shared `../trex_scr_sim_utils.R`. Parts 1–2 print to the console only; Part 3 additionally writes `simulation_results/data/d01_screen_trex_mc_n300_p1000.txt` (summary table) and `.csv` (tidy long format). Results are statistically comparable to the C++ reference, not bit-identical (different RNG streams).

## Interpretation

- The single runs sanity-check that Screen-TRex runs end-to-end and returns a sensible selection at a comfortable SNR of 5 before looking at aggregate MC behavior.
- In the recorded 60-run sweep, Ordinary FDR falls from ~0.47 at SNR 0.01 to ~0.06 for SNR >= 1, with TPR rising from near zero to 1.00 at SNR 5; the Estimated FDR tracks the realized FDR reasonably, sitting somewhat above it at high SNR (~0.09 vs. ~0.06).
- Bootstrap-CI is the more conservative rule throughout: lower FDR (~0.03–0.05 at SNR >= 1) but a clear TPR cost (0.40 at SNR 1 and 0.88 at SNR 5, vs. 0.81 and 1.00 for Ordinary). It is also slower per run due to the bootstrap replicates.

**Last updated**: 2026-07-08
