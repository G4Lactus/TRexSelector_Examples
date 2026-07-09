# Demo 02: Monte Carlo Simulation with Fixed Support

## Purpose

Evaluate classical T-Rex FDR/TPR across an SNR sweep with a **fixed** true support,
benchmarking all base solvers under a shared DGP and FDR target. The foundational
empirical validation demo.

## Data-generating process

`dgp_gauss_snr` — i.i.d. Gaussian design, `y = X beta + eps` calibrated to SNR.
High-dimensional config: `n = 300`, `p = 1000`, `s = 10`, `tFDR = 0.1`. The true
support is drawn **once** before the loop via `make_support_random(s, p, seed = 24)`
and reused across every solver, SNR level, and MC trial (it is not the first 10
features); coefficients are fixed `beta_j = 1`. SNR sweep `{0.1, 0.5, 1.0, 2.0, 5.0}`,
`num_MC = 200` per solver x SNR.

## What it computes

- `run_trex_02` runs the MC over `SOLVERS_DEFAULT` (the 14 shared base solvers:
  TLARS, TLASSO, TENET, TSTAGEWISE, TSTEPWISE, TOMP, TGP, TACGP, TMP, TAFS_rho_0.3,
  TAFS_rho_1.0, TNCGMP_v1, TNCGMP_v0, TOOLS) then over SNR.
- Base control: `K = 20`, `max_dummy_multiplier = 10`, `use_max_T_stop = TRUE`,
  `lloop_strategy = "STANDARD"`, `tloop_stagnation_stop = TRUE`,
  `tloop_max_stagnant_steps = 7`, `opt_threshold = 0.75`. Per-solver params are
  added by `.make_trex_ctrl`.
- Reports averaged FDR and TPR per solver x SNR (`track_L_T = FALSE`, so Avg L / Avg T
  are not computed here).

MC trials run in parallel over a hardcoded `num_cores <- 6L` (set at the top of the
file; edit there to change).

## Running

```bash
Rscript R/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support/demo_trex_02_mc_sim_fixed_support.R
```

Results are written to this demo's own `simulation_results/` folder as a `.txt` table
and a tidy `.csv` (header `solver,metric,snr,value`), stem
`demo_trex_02_mc_sim_fixed_support_results_n300_p1000_stagnation_window_7`. The script
sources `../../support_generators.R`, `../../simulation_utils.R`, and `../trex_sim_utils.R`.

**Last updated**: 2026-07-08
