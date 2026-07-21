# Demo 03: Monte Carlo Simulation with Variable Support (Python)

## Purpose

Investigate T-Rex selector FDR, TPR, average dummy multiplier `L`, and average
stopping time `T` when the **positions** of the active variables are redrawn on
every Monte Carlo trial (support size fixed). Isolates the effect of support
location, complementing Demo 02. Python port of
`cpp/trex_selector_methods/trex/demo_trex_03_mc_sim_variable_support`.

## DGP / Data

- `n = 300`, `p = 1000` (high-dimensional), `s = 10`
- Support **redrawn per trial** (`fixed_support=None`) with fixed cardinality 10
- Fixed coefficients `beta_j = 1` (`rnd_coef=False`)
- SNR grid: `{0.1, 0.2, ..., 2.0}` plus `5.0` (21 values)
- `num_MC = 100` trials per solver x SNR (`_NUM_MC = 100`; deliberate Python
  downscale of the C++ `num_MC = 200`)

## What It Computes

`trx_03_mc_sim_variable_support()` sweeps every solver in `SOLVERS_DEFAULT` (14
solvers) x each SNR via `run_mc_trex`. Base control: `K = 20`,
`max_dummy_multiplier = 10`, `lloop_strategy = "HCONCAT"`,
`tloop_stagnation_stop = True`, `tloop_max_stagnant_steps = 5`. With
`track_L_T=True` it records mean FDR, TPR, Avg L, and Avg T per solver x SNR.

## Imports

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `trex_sim_common` module (explicit `from`-imports of its
symbols); the workers build the selector from `trex_selector_neo` internally —
this demo file itself needs no root alias.

## Output

Writes results to this demo's own `simulation_results/data/` folder via
`save_mc_results`, with stem
`demo_trex_03_mc_sim_variable_support_trex_results_n300_p1000_stagnation_window_5`
(aligned `.txt` table plus tidy `.csv` with FDR/TPR/AvgL/AvgT rows; same stem
convention as the C++ result files).

## Running

```bash
python Python/trex_selector_methods/trex/demo_trex_03_mc_sim_variable_support/demo_trex_03_mc_sim_variable_support.py
```

The MC loop is parallelized with Python `multiprocessing` (spawn start method)
using `_NUM_WORKERS = 6`.

**Last updated**: 2026-07-21
