# Demo 03: Monte Carlo Simulation with Variable Support

## Purpose

Test robustness of classical T-Rex to the **location** of the active set: the support
indices are redrawn on every MC trial while the cardinality stays fixed. Complements
the fixed-support study in Demo 02 and additionally tracks the dummy multiplier `L`
and stopping time `T`.

## Data-generating process

`dgp_gauss_snr` — i.i.d. Gaussian design, `y = X beta + eps` calibrated to SNR.
High-dimensional config: `n = 300`, `p = 1000`, `s = 10`, `tFDR = 0.1`. Each trial
redraws the support via `make_support_random(s, p, seed = trial_seed)`; coefficients
are fixed `beta_j = 1` (`rnd_coef = FALSE`; an isolated `rnorm(s)` stream is available
when `rnd_coef = TRUE`). Fine SNR grid `seq(0.1, 2.0, by = 0.1)` plus `5.0` (21 values),
`num_MC = 200` per solver x SNR.

## What it computes

- `run_trex_03` runs the MC over the same 14 `SOLVERS_DEFAULT` base solvers as Demo 02,
  then over SNR.
- Base control: `K = 20`, `max_dummy_multiplier = 10`, `use_max_T_stop = TRUE`,
  `lloop_strategy = "HCONCAT"`, `tloop_stagnation_stop = TRUE`,
  `tloop_max_stagnant_steps = 5`.
- Reports averaged FDR, TPR, Avg L, and Avg T per solver x SNR (`track_L_T = TRUE`).

MC trials run in parallel over a hardcoded `num_cores <- 6L` (set at the top of the file).

## Running

```bash
Rscript R/trex_selector_methods/trex/demo_trex_03_mc_sim_variable_support/demo_trex_03_mc_sim_variable_support.R
```

Results are written to this demo's own `simulation_results/` folder as a `.txt` table
and a tidy `.csv` (header `solver,metric,snr,value`, with `FDR`/`TPR`/`AvgL`/`AvgT`
rows), stem `demo_trex_03_mc_sim_variable_support_results_n300_p1000_stagnation_window_5`.
The script sources `../../support_generators.R`, `../../simulation_utils.R`, and
`../trex_sim_utils.R`.

**Last updated**: 2026-07-08
