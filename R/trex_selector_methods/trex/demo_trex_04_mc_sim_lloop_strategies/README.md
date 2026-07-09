# Demo 04: L-Loop Strategy Comparison

## Purpose

Compare the T-Rex **L-loop strategies** — the mechanisms by which the dummy-variable
structure is generated/expanded across L-loop iterations — while holding the base
solver fixed at **TLARS**. Shows how the strategy affects FDR, TPR, and the calibrated
dummy multiplier `L` / stopping time `T`.

## Data-generating process

`dgp_gauss_snr` — i.i.d. Gaussian design, `y = X beta + eps` calibrated to SNR.
High-dimensional config: `n = 300`, `p = 1000`, `s = 10`, `tFDR = 0.1`, fixed
coefficients `beta_j = 1`. The support is a random size-`s` subset redrawn per trial
via `make_support_random(s, p, seed = trial_seed)` (the `block_support = TRUE`
contiguous-block variant is available but its call is commented out). SNR grid
`seq(0.1, 2.0, by = 0.1)` plus `5.0` (21 values), `num_MC = 200` per strategy x SNR.

## Strategies compared

`L_STRATEGIES` — **9 rows** (vs the 8 in the C++ demo; the R port adds `SKIPL_50p`):

- STANDARD — fresh i.i.d. dummy matrix each L-loop iteration.
- HCONCAT — horizontally concatenated dummy columns.
- PERMUTATION — re-use the base dummy matrix via permutations.
- PERMUTATION_DIRECT — seed-based permutations, no base matrix in memory.
- DIRECT — seed-based i.i.d. draws, no base matrix in memory.
- SKIPL_5p / SKIPL_10p / SKIPL_20p / SKIPL_50p — skip the L-loop, fixed
  `L = max_dummy_multiplier * p` at 5p / 10p / 20p / 50p.

Shared control per strategy: `solver = "TLARS"`, `K = 20`, `use_max_T_stop = TRUE`,
`lloop_strategy` swept. Reports averaged FDR, TPR, Avg L, Avg T per strategy x SNR
(`track_L_T = TRUE`). MC trials run in parallel over a hardcoded `num_cores <- 6L`.

## Running

```bash
Rscript R/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies/demo_trex_04_mc_sim_lloop_strategies.R
```

Results are written to this demo's own `simulation_results/` folder as a `.txt` table
and a tidy `.csv` (the "solver" column holds the strategy name), stem
`demo_trex_04_lloop_strategies_results_n300_p1000_random_support`. The script sources
`../../support_generators.R`, `../../simulation_utils.R`, and `../trex_sim_utils.R`.

**Last updated**: 2026-07-08
