# Demo 07: T-Rex+DA+BT on Heavy-Tailed Block + White-Noise Data

## Purpose

Combine Demo 05's white-noise dilution with Demo 06's heavy-tailed robustness study:
heavy-tailed AR(1)-Toeplitz blocks embedded in a larger design padded with i.i.d.
Student-t white-noise columns, with total dimension held fixed. It decomposes how
heavy tails and white-noise dilution each affect FDR control and power for the `BT`
method.

## Data-generating process

`dgp_ht_block_white_snr` builds `M` blocks of size `Q` with heavy-tailed (t(nu))
AR(1)-Toeplitz within-block covariance, plus `p_white` appended i.i.d. t(nu)
white-noise columns; total dimension fixed at `p_total = 500`. Predictors X are
always heavy-tailed; response noise is Gaussian or t(nu) by scenario. Active
variables (`s = M`, one per AR block) lie in the AR part only. Base config: `n=150,
p_total=500, M=5, Q=5 (p_ar=25, p_white=475, s=5), amplitude=1.0, rho=0.8, nu=3.0,
SNR=2.0, tFDR=0.2, K=20, num_MC=200, seed=2026`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "BT", ...),
control = trex_control(solver = "TLARS", K = 20))`, driven through the
`.run_mc_ht_block_white()` helper (`trex_method = "trex+DA+BT"`). Outer loops over
noise scenario `{FALSE = Gaussian noise, TRUE = heavy noise}` x linkage
`{"single", "complete", "average"}`.

## What it runs

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}` with `p_ar=25, p_white=475`.
- Part 2: rho sweep `seq(0.0, 0.9, by = 0.1)`.
- Part 3: Q sweep `{5, 10, ..., 50}` (`p_ar = M*Q`, `p_white = p_total - p_ar`).
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}`.
- Part 5: tFDR sweep `seq(0.05, 0.50, by = 0.05)`.

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_07_bt_ht_block_white/demo_trex_da_07_bt_ht_block_white.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_07_bt_ht_block_white/demo_trex_da_07_bt_ht_block_white.R 8    # 8 worker cores
```

Output is redirected to a log file under `simulation_results/` via `sink()`
(comment out the `sink()` lines to see console output). The script sources the
consolidated `../trex_da_dgps.R` and the shared `../../simulation_utils.R` +
`../../support_generators.R`.

## Interpretation

- Unlike Demo 06 (heavy-tailed blocks without white-noise dilution), FDR here tends to decrease as SNR increases and end up quite low, closer to Demo 05's white-noise-diluted pattern — the large white-noise padding helps FDR control even under heavy tails.
- TPR typically reaches higher levels than Demo 06's plateau at the same nominal SNR, again consistent with the white-noise dilution making the true AR(1) blocks comparatively easier to isolate than in Demo 06's smaller all-block design.
- Read as the combined-effects check: compare against Demo 04 (Gaussian, no white noise), Demo 05 (Gaussian, with white noise), and Demo 06 (heavy-tailed, no white noise) to decompose how heavy tails and white-noise dilution each contribute; the tFDR sweep (Part 5) checks target tracking on this diluted design.

**Last updated**: 2026-07-08
