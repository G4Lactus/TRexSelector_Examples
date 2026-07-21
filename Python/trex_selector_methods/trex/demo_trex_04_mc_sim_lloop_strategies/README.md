# Demo 04: L-Loop Strategy Comparison (Python)

## Purpose

Compare the T-Rex **L-loop strategies** (how the dummy-variable structure is
generated/expanded across L-loop iterations) and repeat the whole comparison
across **three base solvers** (TLARS, TOMP, TAFS), showing their effect on FDR,
TPR, calibrated dummy multiplier `L`, and stopping time `T` — and, crucially,
whether the fixed-budget `SKIPL` FDR overshoot seen under TLARS is specific to
the LARS path or persists under the greedy solvers. Python port of
`cpp/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies` (which
extends the R version: R fixes the base solver at TLARS).

## DGP / Data

- `n = 300`, `p = 1000` (high-dimensional), `s = 10`
- Random support redrawn per trial by default (`block_support=False`); a
  contiguous block support `{0, ..., s-1}` is available via `block_support=True`
- Fixed coefficients `beta_j = 1` (`rnd_coef=False`)
- SNR grid: `{0.1, 0.2, ..., 2.0}` plus `5.0` (21 values)
- `num_MC = 25` trials per strategy x SNR (`_NUM_MC = 25`; deliberate Python
  downscale of the C++ `num_MC = 200` — the full grid is 3 solvers x
  8 strategies x 21 SNR values)
- Base solvers (`_SOLVERS`, outer sweep, one result file pair each): TLARS,
  TOMP, TAFS (`rho_afs = 0.3`)

## L-Loop Strategies

`L_STRATEGIES` defines the **eight rows** (mirrors `make_lloop_strategies()`):

1. **STANDARD** — stored; fresh i.i.d. dummy matrices at each L-loop iteration
   (conservative default).
2. **HCONCAT** — stored; horizontally expand (concatenate) dummy columns;
   prefix-stable.
3. **PERMUTATION** — stored base dummy matrix, re-used via deterministic row
   permutations per experiment.
4. **PERMUTATION_ONDEMAND** — seed-derived base + row permutations; nothing
   stored. Bit-identical to PERMUTATION for the same seed.
5. **ONDEMAND** — seed-derived independent dummies; nothing stored.
6. **SKIPL_5p / SKIPL_10p / SKIPL_20p** — skip the L-loop entirely and use the
   fixed `L = max_dummy_multiplier` in `{5, 10, 20} x p`.

(The former `DIRECT` / `PERMUTATION_DIRECT` names were renamed to `ONDEMAND` /
`PERMUTATION_ONDEMAND`; the old names are no longer accepted by the binding.)

## What It Computes

`trx_04_lloop_strategies()` sweeps solver (outer) x strategy x SNR via
`run_mc_trex` with a per-strategy control (`K = 20`, `max_dummy_multiplier`
from the strategy, `lloop_strategy` from the strategy, `use_max_T_stop = True`).
With `track_L_T=True` it records mean FDR, TPR, Avg L, and Avg T per strategy x
SNR — one result container (and file pair) per base solver.

The C++ reference run (num_MC = 200) shows the adaptive strategies holding FDR
at or below target while the fixed-budget `SKIPL_10p`/`SKIPL_20p` levels
overshoot at high SNR under TLARS; the TOMP/TAFS sweeps test whether that
overshoot is LARS-specific.

## Imports

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `trex_sim_common` module (explicit `from`-imports of its
symbols); the workers build the selector from `trex_selector_neo` internally —
this demo file itself needs no root alias.

## Output

Writes results to this demo's own `simulation_results/data/` folder via
`save_mc_results`, one file pair per base solver, stems
`demo_trex_04_lloop_strategies_results_n300_p1000_{tlars,tomp,tafs}_random_support`
(aligned `.txt` table plus tidy `.csv`; the "solver" CSV column holds the
strategy name — the base solver is encoded in the file name, matching the C++
result files).

## Running

```bash
python Python/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies/demo_trex_04_mc_sim_lloop_strategies.py
```

The MC loop is parallelized with Python `multiprocessing` (spawn start method)
using `_NUM_WORKERS = 6`. The full sweep is large — expect a long run even at
the downscaled `num_MC = 25`.

**Last updated**: 2026-07-21
