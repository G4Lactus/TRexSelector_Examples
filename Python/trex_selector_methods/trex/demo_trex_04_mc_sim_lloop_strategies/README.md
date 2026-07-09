# Demo 04: L-Loop Strategy Comparison (Python)

## Purpose

Compare the T-Rex **L-loop strategies** (how the dummy-variable structure is
generated/expanded across L-loop iterations) with the base solver fixed at
TLARS, showing their effect on FDR, TPR, calibrated dummy multiplier `L`, and
stopping time `T`. Python port of
`cpp/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies`.

## DGP / Data

- `n = 300`, `p = 1000` (high-dimensional), `s = 10`
- Random support redrawn per trial by default (`block_support=False`); a
  contiguous block support `{0, ..., s-1}` is available via `block_support=True`
- Fixed coefficients `beta_j = 1` (`rnd_coef=False`)
- SNR grid: `{0.1, 0.2, ..., 2.0}` plus `5.0` (21 values)
- `num_MC = 25` trials per strategy x SNR (`_NUM_MC = 25`)
- Base solver: TLARS (fixed across strategies)

## L-Loop Strategies

`L_STRATEGIES` defines the strategy rows. The full set (STANDARD, HCONCAT,
PERMUTATION, PERMUTATION_DIRECT, DIRECT, and SKIPL at 5p/10p/20p/...) is present
in the file, but as committed only the SKIPL variants **SKIPL_20p** (`L = 20p`)
and **SKIPL_50p** (`L = 50p`) are active; the others are commented out. Uncomment
entries to widen the sweep.

## What It Computes

`trx_04_lloop_strategies()` sweeps each active strategy x SNR via `run_mc_trex`
with a per-strategy control (`K = 20`, `max_dummy_multiplier` from the strategy,
`lloop_strategy` from the strategy). With `track_L_T=True` it records mean FDR,
TPR, Avg L, and Avg T per strategy x SNR.

## Imports

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `trex_sim_common` module. The package imports as
`trex_selector_neo`.

## Output

Writes results to this demo's own `simulation_results/` folder via
`save_mc_results`, with stem
`demo_trex_04_lloop_strategies_results_n300_p1000_random_support` (aligned `.txt`
table plus tidy `.csv`; the "solver" column holds the strategy name).

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies/demo_trex_04_mc_sim_lloop_strategies.py
```

The MC loop is parallelized with Python `multiprocessing` (spawn start method)
using `_NUM_WORKERS = 6`.

**Last updated**: 2026-07-08
