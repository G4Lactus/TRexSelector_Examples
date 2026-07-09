# Demo 05: T-Rex+DA+BT on Block AR(1) + White-Noise Data

## Purpose

Extend Demo 04's block-diagonal AR(1) design by appending a large block of i.i.d.
white-noise columns, testing whether `BT`-based DA-TRex can still isolate the true
AR(1) blocks when they are a small fraction of a much larger, otherwise-uncorrelated
design. Holding the total dimension fixed lets M and Q vary without changing p_total.

## Data-generating process

`dgp_ar1_block_white_snr` builds `M` AR(1) blocks of size `Q` (`p_ar = M * Q`),
plus `p_white` appended i.i.d. N(0,1) columns, with total dimension fixed at
`p_total = p_ar + p_white = 1000`. Active variables (`s = M`, one per AR block) lie
only in the AR part. Base config: `n=300, p_total=1000, M=5, Q=5 (p_ar=25,
p_white=975, s=5), amplitude=1.0, rho=0.7, SNR=2.0, tFDR=0.1, K=20, num_MC=200,
seed=2026`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "BT", ...),
control = trex_control(solver = "TLARS", K = 20))`, driven through the
`.run_mc_ar1_block_white()` helper (`trex_method = "trex+DA+BT"`). Each part loops
over `hc_dist in {"single", "complete", "average"}`.

## What it runs

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}` with `p_ar=25, p_white=975` fixed, outer loop over linkage.
- Part 2: rho sweep `seq(0.0, 0.9, by = 0.1)`, outer loop over linkage.
- Part 3: Q sweep `{5, 10, ..., 50}`; `p_ar = M*Q` varies 25..250, `p_white = 1000 - p_ar`, `s = M = 5`.
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}`; `p_ar`, `p_white`, and `s = M` all vary.

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_05_bt_ar1_block_white/demo_trex_da_05_bt_ar1_block_white.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_05_bt_ar1_block_white/demo_trex_da_05_bt_ar1_block_white.R 8    # 8 worker cores
```

Output is console-only; the script sources the consolidated `../trex_da_dgps.R` and
the shared `../../simulation_utils.R` + `../../support_generators.R`.

## Interpretation

- FDR is typically very tightly controlled here (well below the lower tFDR=0.1 target) and lower than Demo 04's block-only design: the many uncorrelated white-noise columns dilute the design, leaving much less correlated, false-positive-prone structure to latch onto.
- TPR reaches a strong plateau similar to Demo 04, suggesting the added white-noise columns do not meaningfully hurt detection of the true AR(1) blocks, consistent with T-Rex's robustness to null variables.
- Compare directly against Demo 04 (no white noise, same block core) to isolate the effect of a larger, mostly-null p; the Q/M/linkage sweeps remain the key diagnostics for block size, block count, and linkage choice.

**Last updated**: 2026-07-08
