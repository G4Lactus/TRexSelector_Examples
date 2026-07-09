# Demo: T-ENET (Elastic Net) Solver

## Purpose

This demo illustrates the **T-ENET (Terminating Elastic Net)** solver.

T-ENET combines L1-style sparsity with an additional L2 component, making it useful in settings with correlated predictors and grouped signal structure. The demo shows how the solver is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

---

## Key idea

Like the other terminating solvers, T-ENET stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TENET}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}, \lambda_2).
$$

In this demo, the solver is an R6 class from the TRexSelectorNeo package, constructed with a fixed L2 parameter and then advanced through calls of the form

```r
solver <- TENET_Solver$new(X, D, y, lambda2, normalize = TRUE, intercept = TRUE, verbose = TRUE)
solver$execute_step(T_stop, early_stop = TRUE)
```

so, with `early_stop = TRUE`, the run terminates early once the number of dummy variables in the active set reaches $T_{\text{stop}}$.

---

## Demo structure

### 1. Basic T-Elastic Net with early stopping

The first part runs T-ENET with **internal normalization** and an intercept on synthetic data, in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables, selection-quality diagnostics, and a T-ENET-specific diagnostics line with the number of removals and the cycling ratio (`get_num_removals()` / `get_cycling_ratio()`).

### 2. External normalization

The second part applies external L2 normalization via `LpNormScaler` to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-ENET with `normalize = FALSE` and `intercept = FALSE`.

This demonstrates the externally normalized workflow and allows it to be compared conceptually with the internally normalized workflow from the first part. The removals and cycling-ratio diagnostics are printed again.

### 3. Serialization and path consistency

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to `tenet_checkpoint.bin` in the demo folder, reloads it, and continues from the checkpoint to the same final stopping level.

Unlike the cpp version, which uses a static `load(filename, X, D)` factory, the R binding hydrates an existing solver instance: a fresh solver is constructed with the same `X`, `D`, `y`, and `lambda2`, and the checkpoint is loaded via `solver$load(path)`.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, equality of the recorded action path, and the removal and cycling-ratio diagnostics. The checkpoint file is removed afterward.

### 4. Memory-mapped data

The fourth part writes the design matrix $X$ to disk via `convert_to_memory_mapped()` and runs T-ENET on the resulting `mmap_matrix` view. The R binding accepts a memory-mapped $X$ while $D$ and $y$ stay in memory; the cpp demo maps all three.

The cpp counterpart generates $p = 500{,}000$ predictors directly on disk; the R port scales this down to $p = 5{,}000$ (with $n = 1{,}000$ and $10p$ dummies). The file `demo_tenet_X.bin` is written into the demo folder and removed at the end. The removals and cycling-ratio diagnostics are printed here as well.

---

## Running

Run from the repository root:

```bash
Rscript R/tsolvers/demo_ts_04_tenet/demo_ts_04_tenet.R
```

The demo requires the TRexSelectorNeo R package. It prints its diagnostics to the console, including configuration summaries, selected variables, selection-quality measures, removal counts, cycling ratios, serialization checks, and memory-mapped execution status messages.

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

For the dedicated equivalence check between the Gram-based T-ENET formulation and the augmented-LASSO formulation, see [validation_ts_01_tenet_aug_comparison](../validation/validation_ts_01_tenet_aug_comparison/README.md).

---

## Notes

- All parts of this demo use a fixed elastic-net L2 parameter `lambda2 = 0.5` in the constructor.
- Synthetic data come from the suite-level helper `ts_demo_utils.R::gen_synthetic_data()`, an R replication of the cpp `datagen::SyntheticData` generator (SNR-calibrated Gaussian noise, seed 42), so the selected path can be compared against ground truth.
- Selection diagnostics come from `print_selection` / `print_selection_quality`, which compute TPP/FDP/precision/recall/F1 via the TRexSelectorNeo `compute_*` helpers. The active-dummy count is reconstructed from the signed action path because the R binding's `get_actives()` returns predictors only.
- The serialization example writes a temporary checkpoint file named `tenet_checkpoint.bin` and removes it afterward.
- The memory-mapped example writes a temporary file for $X$ only, runs T-ENET on the mapped view, and removes the file at the end.

---

**Last updated**: 2026-07-08
