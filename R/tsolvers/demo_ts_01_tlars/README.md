# Demo: T-LARS Solver

## Purpose

This demo illustrates the **T-LARS (Terminating Least Angle Regression)** solver, one of the core terminating solvers used in the T-Rex framework.

It shows how T-LARS is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

The demo is intended to explain solver behavior and interface usage in isolation, outside the full T-Rex calibration loops.

---

## Key idea

A terminating solver does not compute a full path over all possible steps. Instead, it stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TLARS}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, the solver is an R6 class from the TRexSelectorNeo package, and early stopping is triggered through calls of the form:

```r
solver <- TLARS_Solver$new(X, D, y, normalize = TRUE, intercept = TRUE, verbose = TRUE)
solver$execute_step(T_stop, early_stop = TRUE)
```

where `T_stop` is the number of dummy variables that must enter the active set before the solver terminates early: with `early_stop = TRUE` the run stops as soon as the count of active dummies reaches $T_{\text{stop}}$.

---

## Demo structure

### 1. Basic T-LARS with early stopping

The first part runs T-LARS with **internal normalization** and an intercept on synthetic data, and it is executed in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support (1-based indices $\{28, 150, 399, 421, 5\}$).

### 2. External normalization

The second part preprocesses the signal matrix and dummy matrix externally using `LpNormScaler` (L2-based column normalization with centering), centers the response manually, and then runs T-LARS with `normalize = FALSE` and `intercept = FALSE`.

This section is useful for comparing the solver's external-normalization workflow against the internally normalized workflow from the first demo.

### 3. Serialization and warm start

The third part runs a reference T-LARS fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to `tlars_checkpoint.bin` in the demo folder, reloads it, and continues to the same final stopping level.

Unlike the cpp version, which uses a static `load(filename, X, D)` factory, the R binding hydrates an existing solver instance: a fresh solver is constructed with the same `X`, `D`, and `y`, and the checkpoint is loaded via `solver$load(path)`.

The reloaded solver is compared against the reference run through the number of steps, the final RSS difference, the final $R^2$ difference, and the equality of the recorded action path. The checkpoint file is removed at the end.

### 4. Memory-mapped data

The fourth part writes the design matrix $X$ to disk via `convert_to_memory_mapped()` and runs T-LARS on the resulting `mmap_matrix` view. The R binding accepts a memory-mapped $X$ while $D$ and $y$ stay in memory; the cpp demo maps all three.

The cpp counterpart generates $p = 500{,}000$ predictors directly on disk; the R port scales this down to $p = 5{,}000$ (with $n = 1{,}000$ and $10p$ dummies). The file `demo_tlars_X.bin` is written into the demo folder and removed at the end.

---

## Running

Run from the repository root:

```bash
Rscript R/tsolvers/demo_ts_01_tlars/demo_ts_01_tlars.R
```

The demo requires the TRexSelectorNeo R package. It prints its diagnostics to the console, including configuration summaries, selected variables, selection-quality measures, serialization checks, and memory-mapped execution status messages.

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

- Synthetic data come from the suite-level helper `ts_demo_utils.R::gen_synthetic_data()`, an R replication of the cpp `datagen::SyntheticData` generator (SNR-calibrated Gaussian noise, seed 42), so the selected path can be compared against ground truth.
- Selection diagnostics come from `print_selection` / `print_selection_quality`, which compute TPP/FDP/precision/recall/F1 via the TRexSelectorNeo `compute_*` helpers. The active-dummy count is reconstructed from the signed action path because the R binding's `get_actives()` returns predictors only.
- The serialization example writes a temporary checkpoint file and removes it after the comparison.
- The memory-mapped example writes a temporary file for $X$ only, runs the solver on the mapped view, and then removes the file at the end.

---

**Last updated**: 2026-07-08
