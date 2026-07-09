# Demo 02: Monte Carlo Simulation with Fixed Support (Python)

## Purpose

Evaluate T-Rex selector FDR and TPR across a range of Signal-to-Noise Ratios
with a fixed true support, comparing all default base solvers. Python port of
`cpp/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support`.

## DGP / Data

- `n = 300`, `p = 1000` (high-dimensional), `s = 10`
- True support drawn **once** via `make_support_random(s, p, seed=24)` and shared
  across every solver, SNR level, and MC trial
- Fixed coefficients `beta_j = 1`
- SNR sweep: `{0.1, 0.5, 1.0, 2.0, 5.0}`
- `num_MC = 100` trials per solver x SNR (`_NUM_MC = 100`)

## What It Computes

`trx_02_mc_sim_fixed_support()` sweeps every solver in `SOLVERS_DEFAULT` (the
14-solver default list) x each SNR, calling `run_mc_trex` per cell. Base control:
`K = 20`, `max_dummy_multiplier = 10`, `lloop_strategy = "STANDARD"`,
`tloop_stagnation_stop = True`, `tloop_max_stagnant_steps = 7`,
`opt_threshold = 0.75`, `use_memory_mapping = False`. It records mean FDR and TPR
per solver x SNR (`track_L_T=False`).

## Imports

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `support_generators` and `trex_sim_common` modules. The
package imports as `trex_selector_neo`.

## Output

Writes results to this demo's own `simulation_results/` folder via
`save_mc_results`, with stem
`demo_trex_02_mc_sim_fixed_support_results_n300_p1000_stagnation_window_7`
(aligned `.txt` table plus tidy `.csv`).

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support/demo_trex_02_mc_sim_fixed_support.py
```

The MC loop is parallelized with Python `multiprocessing` (spawn start method)
using `_NUM_WORKERS = 6`.

**Last updated**: 2026-07-08
