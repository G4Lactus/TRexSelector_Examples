# Demo: T-LASSO Solver (Python)

## Purpose

This demo illustrates the **T-LASSO (Terminating Least Absolute Shrinkage and Selection Operator)** solver.

It shows how T-LASSO is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

In addition to the standard selection diagnostics, the demo also reports T-LASSO-specific path diagnostics such as variable removals and cycling behavior. It mirrors `cpp/tsolvers/demo_ts_02_tlasso`.

---

## Key idea

Like the other terminating solvers, T-LASSO stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set, rather than tracing a full path:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TLASSO}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form

```python
solver.executeStep(T_stop, early_stop=True)
```

where `T_stop` is the number of dummy variables that must enter the active set before the solver terminates early: with `early_stop=True` the run stops as soon as the count of active dummies reaches $T_{\text{stop}}$.

The solver class is `TLASSO_Solver` from `trex_selector_neo.tsolvers.lars_based`.

---

## Demo structure

### 1. Basic T-LASSO with early stopping

The first part runs T-LASSO with **internal normalization** and an intercept on synthetic data. It is executed in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables, selection-quality diagnostics, and T-LASSO-specific quantities: the number of removals (`getNumRemovals()`) and the cycling ratio (`getCyclingRatio()`).

### 2. External normalization

The second part applies external L2 normalization (`LpNormScaler` with `NormType.L2`) to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-LASSO with internal normalization disabled.

This section demonstrates the externally normalized workflow and again reports the resulting selection path, selection quality, removals, and cycling statistics.

### 3. Serialization and path consistency

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues from the checkpoint to the same final stopping level.

Unlike cpp, which uses a static `load(filename, X, D)` factory, the Python binding hydrates an existing solver instance: a new `TLASSO_Solver` is constructed with the same `X`, `D`, and `y`, and then `solver.load(path)` restores the checkpointed state.

The reloaded path is compared against the reference run using the final number of steps, the final RSS difference, the final $R^2$ difference, equality of the recorded action path, and the T-LASSO-specific removal and cycling diagnostics.

### 4. Memory-mapped data

The fourth part writes $X$, $D$, and $y$ to disk with `numpy_to_memmap()` and runs T-LASSO on the zero-copy mmap views (via `.to_numpy()`), again reporting the removal and cycling diagnostics.

Note: the cpp demo generates $p = 500{,}000$ on disk; the Python port scales this down to $p = 5{,}000$ (with $n = 1000$ and $10p$ dummies).

---

## Running

The demo is a script — run it from the repo root with:

```bash
.venv/bin/python Python/tsolvers/demo_ts_02_tlasso/demo_ts_02_tlasso.py
```

There is no build step; it requires the `trex_selector_neo` package installed in the repo `.venv`.

The demo prints its diagnostics to the console, including configuration summaries, selected variables, path-quality measures, removal counts, cycling ratios, serialization checks, and memory-mapped execution status messages.

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

---

## Notes

- Synthetic data come from the suite-level helper `ts_demo_utils.gen_synthetic_data()` (numpy replication of the cpp `datagen::SyntheticData`, SNR-calibrated Gaussian noise, seed 42), so the selection path can be compared against ground truth; diagnostics use `print_selection` / `print_selection_quality` (TPP/FDP/precision/recall/F1 via `trex_selector_neo.utils`), with 0-based support indices matching cpp.
- The removals/cycling-ratio diagnostics block is printed after each of the four parts.
- The thread count is set to 6 via `trex_selector_neo.utils.set_num_threads`, mirroring the cpp OpenMP setup.
- The serialization example writes a temporary checkpoint file `tlasso_checkpoint.bin` into the demo folder and removes it afterward.
- The memory-mapped example writes `demo_tlasso_{X,D,y}.bin` into the demo folder, runs T-LASSO on those mapped views, and removes those files at the end.
- For a dedicated cross-check against related LARS/LASSO reference behavior, see the validation folder under `Python/tsolvers/validation`.

---

**Last updated**: 2026-07-08
