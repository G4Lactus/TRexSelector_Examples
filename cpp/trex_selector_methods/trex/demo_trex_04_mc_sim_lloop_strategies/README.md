# Demo 04: L-Loop Strategy Comparison

## Purpose

Compare the **L-loop strategies** of the T-Rex Selector — the different mechanisms by which the dummy-variable structure is generated/expanded across L-loop iterations — while holding the base solver fixed at **TLARS**. This shows how the L-loop strategy affects FDR control, power (TPR), and the calibrated dummy multiplier $L$ / stopping time $T$.

---

## Data Generation Parameters

- **Sample size**: $n = 300$
- **Number of features**: $p = 1000$
- **True support**: a random subset of size $s = 10$, **redrawn per trial** (`block_support = false`) by shuffling $\{0, \ldots, p-1\}$ with `std::mt19937(seed + 500000)`, keeping 10 indices and sorting them. A contiguous block support $\{0, \ldots, s-1\}$ is available via the `block_support = true` call, which is compiled but guarded by `if (false)`.
- **True coefficients**: fixed $\beta_j = 1$ (`rnd_coef = false`)
- **SNR grid**: 21 values — $\{0.1, 0.2, \ldots, 2.0\}$ plus $5.0$
- **Monte Carlo repetitions**: `num_MC = 10` trials per strategy × SNR level
- **DGP**: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, \sigma^2 I_n)$, Normal predictors and Normal noise
- **Base solver**: TLARS (fixed across all strategies)
- **tFDR**: $0.1$

---

## L-Loop Strategies Tested

`make_lloop_strategies()` sweeps the **six L-loop strategy types** defined by the library `LLoopStrategy` enum, with `SKIPL` evaluated at three fixed dummy levels (5p, 10p, 20p) — **eight rows** in total.

The strategies span two orthogonal axes: **dummy source** (fresh independent draws per experiment vs one shared base matrix whose rows are permuted per experiment) and **storage** (*stored* strategies keep matrices in the `DummyGenerator` for the whole run; *on-demand* strategies re-derive everything from the seed at each step, with zero persistent state):

1. **STANDARD** — stored; fresh i.i.d. dummy matrices at each L-loop iteration (conservative default, matches the CRAN R reference).
2. **HCONCAT** — stored; horizontally expand (concatenate) dummy columns; prefix-stable.
3. **PERMUTATION** — stored base dummy matrix, re-used via deterministic **row permutations** per experiment.
4. **PERMUTATION_ONDEMAND** — seed-derived base + row permutations per experiment; nothing stored. Bit-identical to PERMUTATION for the same seed.
5. **ONDEMAND** — seed-derived independent dummies per experiment; nothing stored.
6. **SKIPL** — skip the L-loop entirely and use $L = \text{max\_dummy\_multiplier}$, evaluated at:
   - **SKIPL_5p** ($L = 5p$)
   - **SKIPL_10p** ($L = 10p$)
   - **SKIPL_20p** ($L = 20p$)

The library enum has **no** "Doubling" or "Fixed" strategy. (The former `DIRECT` / `PERMUTATION_DIRECT` names were renamed to `ONDEMAND` / `PERMUTATION_ONDEMAND`; the old names no longer exist.)

---

## Control Parameters

Set per strategy (see `make_lloop_strategies()` / the strategy loop):

```
solver_type = TLARS              # Base solver, fixed
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Default for adaptive strategies; 5 / 10 / 20 for SKIPL
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
dummy_distribution = Normal      # Dummy predictors drawn from N(0,1)
lloop_strategy = <varies>        # Main variable swept across the 6 types
tFDR = 0.1                       # Target FDR control level
```

The MC loop is parallelized with OpenMP (`omp_set_num_threads(6)`).

---

## Output Files

Both files are written to `simulation_results/data/`. The stem encodes the support scenario (`random_support` for the active run):

### Main Result File
**`demo_trex_04_lloop_strategies_results_n300_p1000_random_support.txt`**

Aligned table (written by the shared `save_and_print_mc_results`) with four metric rows — FDR, TPR, Avg L, Avg T — per strategy across the 21 SNR columns. The row label is the strategy name (STANDARD, HCONCAT, PERMUTATION, PERMUTATION_ONDEMAND, ONDEMAND, SKIPL_5p, SKIPL_10p, SKIPL_20p):

```
======================================================================
=== T-Rex Results (averaged over 10 Monte Carlo runs) ===
======================================================================

Solver         Metric    SNR       0.1       0.2  ...       2.0       5.0
--------------------------------------------------------------------------
STANDARD       FDR              ...
               TPR              ...
               Avg L            ...
               Avg T            ...

HCONCAT        ...
...
```

### Tidy-Format CSV
**`demo_trex_04_lloop_strategies_results_n300_p1000_random_support.csv`**

Long/stacked format, header column order **`solver,metric,snr,value`** (here the "solver" column holds the strategy name), with `FDR`, `TPR`, `AvgL`, `AvgT` rows.

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies/demo_trex_04_mc_sim_lloop_strategies
```

---

## Key Questions Addressed

1. **Do all L-loop strategies maintain FDR control?**
   - Expected: FDR $\leq$ tFDR across strategies.

2. **How does the calibrated $L$ differ between adaptive strategies and the SKIPL levels?**
   - The Avg L rows show the adaptive strategies' calibrated $L$ versus SKIPL's fixed 5p/10p/20p.

3. **Is TPR comparable across strategies at matched SNR?**
   - Expected: broadly comparable power, with strategy-specific variation.

---

## Interpretation Guide

**What to look for:**
- **FDR equivalence**: all strategies should keep FDR $\leq$ tFDR.
- **TPR consistency**: minor variation due to randomness and the small MC count (`num_MC = 10`).
- **Avg L / Avg T**: how each strategy calibrates the dummy structure and when the T-loop stops.

**Practical significance:**
- The seed-based ONDEMAND / PERMUTATION_ONDEMAND strategies avoid holding any dummy matrix in memory — the preferred choice at scale.
- PERMUTATION and PERMUTATION_ONDEMAND produce bit-identical selections for the same seed; they differ only in memory footprint.
- SKIPL trades adaptivity for a fixed, larger dummy budget.

This C++ demo mirrors `R/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies.R`.

---

**Last updated**: 2026-07-10
