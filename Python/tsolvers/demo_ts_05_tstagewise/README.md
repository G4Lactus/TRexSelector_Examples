# Demo: T-Stagewise Solver (Python)

## Purpose

This demo illustrates the **T-Stagewise (Terminating Incremental Forward Stagewise Regression)** solver.

T-Stagewise proceeds through small incremental updates in directions chosen from the current residual correlations, which typically produces a more gradual path than more aggressive greedy methods. The demo shows how the solver is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

It mirrors `cpp/tsolvers/demo_ts_05_tstagewise`.

---

## Key idea

Like the other terminating solvers, T-Stagewise stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TStagewise}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form

```python
solver.executeStep(T_stop, early_stop=True)
```

so, with `early_stop=True`, the solver terminates early once the number of dummy variables in the active set reaches $T_{\text{stop}}$.

The solver class is `TSTAGEWISE_Solver` from `trex_selector_neo.tsolvers.lars_based`.

---

## Demo structure

### 1. Basic T-Stagewise with early stopping

The first part runs T-Stagewise with **internal normalization** and an intercept on synthetic data, in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables, selection-quality diagnostics, and T-Stagewise-specific path diagnostics: the number of removals (`getNumRemovals()`) and the cycling ratio (`getCyclingRatio()`).

### 2. External normalization

The second part applies external L2 normalization (`LpNormScaler` with `NormType.L2`) to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-Stagewise with internal normalization disabled, again reporting removals and cycling statistics.

### 3. Serialization and warm start

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues from the checkpoint to the same final stopping level.

Unlike cpp, which uses a static `load(filename, X, D)` factory, the Python binding hydrates an existing solver instance: a new `TSTAGEWISE_Solver` is constructed with the same `X`, `D`, and `y`, and then `solver.load(path)` restores the checkpointed state.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, and equality of the recorded action path.

### 4. Memory-mapped data

The fourth part writes $X$, $D$, and $y$ to disk with `numpy_to_memmap()` and runs T-Stagewise on the zero-copy mmap views (via `.to_numpy()`), again reporting the removal and cycling diagnostics.

Note: the cpp demo generates $p = 500{,}000$ on disk; the Python port scales this down to $p = 5{,}000$ (with $n = 1000$ and $10p$ dummies).

---

## Running

The demo is a script — run it from the repo root with:

```bash
.venv/bin/python Python/tsolvers/demo_ts_05_tstagewise/demo_ts_05_tstagewise.py
```

There is no build step; it requires the `trex_selector_neo` package installed in the repo `.venv`.

The demo prints its diagnostics to the console, including configuration summaries, selected variables, selection-quality measures, removal counts, cycling ratios, serialization checks, and memory-mapped execution status messages.

---

## What to look for

When reading the console output, focus on the following points:

- whether the selected variables recover the known true support reasonably well,
- how the solver behaves under both high-dimensional and low-dimensional settings,
- whether externally normalized and internally normalized workflows both execute successfully,
- whether the serialized solver state can be reloaded and continued without changing the final path,
- how many active-set removals occur,
- whether the reported cycling ratio stays well behaved,
- whether the memory-mapped workflow completes successfully and cleans up its temporary files.

Because T-Stagewise uses incremental updates, it may require more steps to build up the same effective support size than more aggressive path-following solvers, and its path is typically smoother than a purely greedy stepwise selection path.

---

## Notes

- Synthetic data come from the suite-level helper `ts_demo_utils.gen_synthetic_data()` (numpy replication of the cpp `datagen::SyntheticData`, SNR-calibrated Gaussian noise, seed 42), so the selected path can be compared against ground truth; selection diagnostics use `print_selection` / `print_selection_quality` (TPP/FDP/precision/recall/F1 via `trex_selector_neo.utils`), with 0-based support indices matching cpp.
- The removals/cycling-ratio diagnostics block is printed after Parts 1, 2, and 4 (the serialization part reports only the path-consistency comparison).
- The thread count is set to 6 via `trex_selector_neo.utils.set_num_threads`, mirroring the cpp OpenMP setup.
- The serialization example writes a temporary checkpoint file `tstagewise_checkpoint.bin` into the demo folder and removes it afterward.
- The memory-mapped example writes `demo_tstagewise_{X,D,y}.bin` into the demo folder, runs T-Stagewise on those mapped views, and removes those files at the end.

---

**Last updated**: 2026-07-08
