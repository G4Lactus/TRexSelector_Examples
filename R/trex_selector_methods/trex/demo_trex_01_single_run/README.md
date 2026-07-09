# Demo 01: Single-Run Classical T-Rex Selector

## Purpose

Demonstrate basic T-Rex selector usage and its output format on a single random
draw, using the TLARS base solver. The demo runs two contrasting scenarios to show
the reporting of the internal T-Rex diagnostics.

## Data-generating process

`dgp_gauss_snr` builds i.i.d. Gaussian design `X ~ N(0,1)` and `y = X beta + eps`
with noise calibrated to a target SNR. Fixed true support
`TRUE_SUPPORT_P1 = c(28, 150, 44, 129, 43, 5)` (1-based; the C++ 0-based set is
`{27, 149, 43, 128, 42, 4}`), fixed coefficients `beta_j = 1` (`amplitude = 1`),
`SNR = 1.0`, `tFDR = 0.1`, `seed = 4231`. Two dimension settings:

- Low-dimensional (`n > p`): `n = 5000`, `p = 1000` — run first.
- High-dimensional (`p > n`): `n = 150`, `p = 300` — run second.

## Control

`trex_control(solver = "TLARS", K = 20, max_dummy_multiplier = 10,
use_max_T_stop = TRUE, lloop_strategy = "HCONCAT", tloop_stagnation_stop = FALSE,
tloop_max_stagnant_steps = 5, use_openmp = FALSE)`.

## What it computes

- Part 1 (`run_part_01`): one `TRexSelector$new(...)$select()` per dimension setting.
- `.print_single_run` reports the selected indices, TPP, FDP, and the internal
  T-Rex statistics `phi_prime`, `phi_mat`, `fdp_hat_mat`, `r_mat`, and the voting grid.

This is not a Monte Carlo study; results are from a single seed.

## Running

```bash
Rscript R/trex_selector_methods/trex/demo_trex_01_single_run/demo_trex_01_single_run.R
```

Output is console-only (no result files are written). The script sources the shared
`../../support_generators.R` and `../trex_sim_utils.R`.

**Last updated**: 2026-07-08
