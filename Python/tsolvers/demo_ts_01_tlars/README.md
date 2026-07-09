# Demo: T-LARS Solver (Python)

## Purpose

This demo illustrates the **T-LARS (Terminating Least Angle Regression)** solver, one of the core terminating solvers used in the T-Rex framework.

It shows how T-LARS is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

The demo is intended to explain solver behavior and interface usage in isolation, outside the full T-Rex calibration loops. It mirrors `cpp/tsolvers/demo_ts_01_tlars`.

---

## Key idea

A terminating solver does not compute a full path over all possible steps. Instead, it stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TLARS}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form:

```python
solver.executeStep(T_stop, early_stop=True)
```

where `T_stop` is the number of dummy variables that must enter the active set before the solver terminates early: with `early_stop=True` the run stops as soon as the count of active dummies reaches $T_{\text{stop}}$.

The solver class is `TLARS_Solver` from `trex_selector_neo.tsolvers.lars_based`.

---

## Demo structure

### 1. Basic T-LARS with early stopping

The first part runs T-LARS with **internal normalization** and an intercept on synthetic data, and it is executed in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part preprocesses the signal matrix and dummy matrix externally using L2-based column normalization (`LpNormScaler` with `NormType.L2`), centers the response manually, and then runs T-LARS with internal normalization disabled.

This section is useful for comparing the solver's external-normalization workflow against the internally normalized workflow from the first demo.

### 3. Serialization and warm start

The third part runs a reference T-LARS fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues to the same final stopping level.

Unlike cpp, which uses a static `load(filename, X, D)` factory, the Python binding hydrates an existing solver instance: a new `TLARS_Solver` is constructed with the same `X`, `D`, and `y`, and then `solver.load(path)` restores the checkpointed state.

The reloaded solver is compared against the reference run through the number of steps, the final RSS difference, the final $R^2$ difference, and the equality of the recorded action path (`getActions()`).

### 4. Memory-mapped data

The fourth part writes $X$, $D$, and $y$ to disk with `numpy_to_memmap()` and runs T-LARS on the zero-copy mmap views (via `.to_numpy()`).

Note: the cpp demo generates $p = 500{,}000$ on disk; the Python port scales this down to $p = 5{,}000$ (with $n = 1000$ and $10p$ dummies).

This demonstrates that the solver can also be used in workflows where the data are too large or too inconvenient to keep entirely in ordinary in-memory arrays.

---

## Running

The demo is a script — run it from the repo root with:

```bash
.venv/bin/python Python/tsolvers/demo_ts_01_tlars/demo_ts_01_tlars.py
```

There is no build step; it requires the `trex_selector_neo` package installed in the repo `.venv`.

The demo prints its diagnostics to the console, including configuration summaries, selected variables, selection-quality measures, serialization checks, and memory-mapped execution status messages.

---

## What to look for

When reading the console output, focus on the following points:

- whether the selected variables recover the known true support reasonably well,
- how the solver behaves under both high-dimensional and low-dimensional settings,
- whether externally normalized and internally normalized workflows both execute successfully,
- whether the serialized solver state can be reloaded and continued without changing the final path,
- whether the memory-mapped workflow completes successfully and cleans up its generated files.

---

## Notes

- Synthetic data come from the suite-level helper `ts_demo_utils.gen_synthetic_data()`, a numpy replication of the cpp `datagen::SyntheticData` (SNR-calibrated Gaussian noise, seed 42), so the selected path can be compared against ground truth.
- Selection diagnostics are printed with `print_selection` / `print_selection_quality` (TPP/FDP/precision/recall/F1 via `trex_selector_neo.utils`). Support indices are 0-based, matching cpp.
- The thread count is set to 6 via `trex_selector_neo.utils.set_num_threads`, mirroring the cpp OpenMP setup.
- The serialization example writes a temporary checkpoint file `tlars_checkpoint.bin` into the demo folder and removes it after the comparison.
- The memory-mapped example writes `demo_tlars_{X,D,y}.bin` into the demo folder, runs the solver on mapped views, and then removes those files at the end.

---

**Last updated**: 2026-07-08
