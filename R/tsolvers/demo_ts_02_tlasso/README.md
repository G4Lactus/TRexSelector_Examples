# Demo: T-LASSO Solver

## Purpose

This demo illustrates the **T-LASSO (Terminating Least Absolute Shrinkage and Selection Operator)** solver.

It shows how T-LASSO is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

In addition to the standard selection diagnostics, the demo also reports T-LASSO-specific path diagnostics, namely the number of variable removals along the path.

---

## Key idea

Like the other terminating solvers, T-LASSO stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set, rather than tracing a full path:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TLASSO}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, the solver is an R6 class from the TRexSelectorNeo package, and early stopping is triggered through calls of the form

```r
solver <- TLASSO_Solver$new(X, D, y, normalize = TRUE, intercept = TRUE, verbose = TRUE)
solver$execute_step(T_stop, early_stop = TRUE)
```

where `T_stop` is the number of dummy variables that must enter the active set before the solver terminates early: with `early_stop = TRUE` the run stops as soon as the count of active dummies reaches $T_{\text{stop}}$.

---

## Demo structure

### 1. Basic T-LASSO with early stopping

The first part runs T-LASSO with **internal normalization** and an intercept on synthetic data. It is executed in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables, selection-quality diagnostics, and a T-LASSO-specific diagnostics line with the number of removals (`solver$get_num_removals()`). Unlike the cpp demo, no cycling ratio is reported, because the R binding does not expose the cycling ratio for T-LASSO.

### 2. External normalization

The second part applies external L2 normalization via `LpNormScaler` to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-LASSO with `normalize = FALSE` and `intercept = FALSE`.

This section demonstrates the externally normalized workflow and again reports the resulting selection path, selection quality, and removal count.

### 3. Serialization and path consistency

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to `tlasso_checkpoint.bin` in the demo folder, reloads it, and continues from the checkpoint to the same final stopping level.

Unlike the cpp version, which uses a static `load(filename, X, D)` factory, the R binding hydrates an existing solver instance: a fresh solver is constructed with the same `X`, `D`, and `y`, and the checkpoint is loaded via `solver$load(path)`.

The reloaded path is compared against the reference run using the final number of steps, the final RSS difference, the final $R^2$ difference, equality of the recorded action path, and the removal count. The checkpoint file is removed afterward.

### 4. Memory-mapped data

The fourth part writes the design matrix $X$ to disk via `convert_to_memory_mapped()` and runs T-LASSO on the resulting `mmap_matrix` view. The R binding accepts a memory-mapped $X$ while $D$ and $y$ stay in memory; the cpp demo maps all three.

The cpp counterpart generates $p = 500{,}000$ predictors directly on disk; the R port scales this down to $p = 5{,}000$ (with $n = 1{,}000$ and $10p$ dummies). The file `demo_tlasso_X.bin` is written into the demo folder and removed at the end.

---

## Running

Run from the repository root:

```bash
Rscript R/tsolvers/demo_ts_02_tlasso/demo_ts_02_tlasso.R
```

The demo requires the TRexSelectorNeo R package. It prints its diagnostics to the console, including configuration summaries, selected variables, path-quality measures, removal counts, serialization checks, and memory-mapped execution status messages.

---

## What to look for

When reading the console output, focus on the following points:

- whether the selected variables recover the known true support reasonably well,
- how the solver behaves under both high-dimensional and low-dimensional settings,
- whether externally normalized and internally normalized workflows both execute successfully,
- whether the serialized solver state can be reloaded and continued without changing the final path,
- how many active-set removals occur,
- whether the memory-mapped workflow completes successfully and cleans up its temporary files.

---

## Notes

- Synthetic data come from the suite-level helper `ts_demo_utils.R::gen_synthetic_data()`, an R replication of the cpp `datagen::SyntheticData` generator (SNR-calibrated Gaussian noise, seed 42), so the selection path can be compared against ground truth.
- Selection diagnostics come from `print_selection` / `print_selection_quality`, which compute TPP/FDP/precision/recall/F1 via the TRexSelectorNeo `compute_*` helpers. The active-dummy count is reconstructed from the signed action path because the R binding's `get_actives()` returns predictors only.
- The demo uses both in-memory and memory-mapped workflows.
- The serialization example writes a temporary checkpoint file named `tlasso_checkpoint.bin` and removes it afterward.
- The memory-mapped example writes a temporary file for $X$ only, runs T-LASSO on the mapped view, and removes the file at the end.
- For a dedicated cross-check against related LARS/LASSO reference behavior, see the validation folder described in `validation_ts_02_tlars_tlasso_rcompare`.

---

**Last updated**: 2026-07-08
