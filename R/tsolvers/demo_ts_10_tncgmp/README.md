# Demo: T-NCGMP (Norm-Corrective Generalized Matching Pursuit) Solver

## Purpose

This demo illustrates the **T-NCGMP (Terminating Norm-Corrective Generalized Matching Pursuit)** solver.

T-NCGMP is a pursuit-style greedy method that augments matching-pursuit-style updates with a norm-correction mechanism, aiming for better numerical stability and improved convergence behavior over longer paths. The demo shows how the solver is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

---

## Key idea

Like the other terminating solvers, T-NCGMP stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TNCGMP}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, the solver is an R6 class from the TRexSelectorNeo package, constructed with the `"line_search"` variant and advanced through calls of the form

```r
solver <- TNCGMP_Solver$new(X, D, y, variant = "line_search", normalize = TRUE,
                            intercept = TRUE, verbose = TRUE)
solver$execute_step(T_stop, early_stop = TRUE)
```

so, with `early_stop = TRUE`, the solver terminates early once the number of dummy variables in the active set reaches $T_{\text{stop}}$.

---

## Demo structure

### 1. Basic T-NCGMP with early stopping

The first part runs T-NCGMP with **internal normalization** and an intercept on synthetic data, in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part applies external L2 normalization via `LpNormScaler` to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-NCGMP with the same solver variant on the transformed data (`normalize = FALSE`, `intercept = FALSE`).

This demonstrates the externally normalized workflow and allows it to be compared conceptually with the internally normalized workflow from the first part.

### 3. Serialization and warm start

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to `tncgmp_checkpoint.bin` in the demo folder, reloads it, and continues from the checkpoint to the same final stopping level. In this part the solvers are constructed with `verbose = TRUE`, so per-step progress is visible in the console.

Unlike the cpp version, which uses a static `load(filename, X, D)` factory, the R binding hydrates an existing solver instance: a fresh solver is constructed with the same `X`, `D`, `y`, and variant, and the checkpoint is loaded via `solver$load(path)`.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, and equality of the recorded action path. The checkpoint file is removed afterward.

### 4. Memory-mapped data

The fourth part writes the design matrix $X$ to disk via `convert_to_memory_mapped()` and runs T-NCGMP on the resulting `mmap_matrix` view. The R binding accepts a memory-mapped $X$ while $D$ and $y$ stay in memory; the cpp demo maps all three.

The cpp counterpart generates $p = 500{,}000$ predictors directly on disk; the R port scales this down to $p = 5{,}000$ (with $n = 1{,}000$ and $10p$ dummies). The file `demo_tncgmp_X.bin` is written into the demo folder and removed at the end.

---

## Running

Run from the repository root:

```bash
Rscript R/tsolvers/demo_ts_10_tncgmp/demo_ts_10_tncgmp.R
```

The demo requires the TRexSelectorNeo R package. It prints its diagnostics to the console, including configuration summaries, selected variables, selection-quality measures, serialization checks, and memory-mapped execution status messages.

---

## What to look for

When reading the console output, focus on the following points:

- whether the selected variables recover the known true support reasonably well,
- how the solver behaves under both high-dimensional and low-dimensional settings,
- whether externally normalized and internally normalized workflows both execute successfully,
- whether the serialized solver state can be reloaded and continued without changing the final path,
- whether the memory-mapped workflow completes successfully and cleans up its temporary files.

Conceptually, T-NCGMP extends the simpler pursuit variants with a norm-corrective step, so it is a natural place to watch for improved stability and more controlled long-path behavior relative to plain T-MP.

---

## Notes

- All parts of this demo construct the solver with `variant = "line_search"`.
- Synthetic data come from the suite-level helper `ts_demo_utils.R::gen_synthetic_data()`, an R replication of the cpp `datagen::SyntheticData` generator (SNR-calibrated Gaussian noise, seed 42), so the selected path can be compared against ground truth.
- Selection diagnostics come from `print_selection` / `print_selection_quality`, which compute TPP/FDP/precision/recall/F1 via the TRexSelectorNeo `compute_*` helpers. The active-dummy count is reconstructed from the signed action path because the R binding's `get_actives()` returns predictors only.
- The serialization example writes a temporary checkpoint file named `tncgmp_checkpoint.bin` and removes it afterward.
- The memory-mapped example writes a temporary file for $X$ only, runs T-NCGMP on the mapped view, and removes the file at the end.

---

**Last updated**: 2026-07-08
