# Demo: T-OLS (Orthogonal Least Squares) Solver

## Purpose

This demo illustrates the **T-OLS (Terminating Orthogonal Least Squares)** solver.

T-OLS is a greedy variable-selection method that, at each step, chooses the variable that gives the largest reduction in residual sum of squares under a least-squares update, rather than selecting purely by raw correlation. The demo shows how the solver is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

---

## Key idea

Like the other terminating solvers, T-OLS stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TOOLS}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, the solver is an R6 class from the TRexSelectorNeo package, and early stopping is triggered through calls of the form

```r
solver <- TOOLS_Solver$new(X, D, y, normalize = TRUE, intercept = TRUE, verbose = TRUE)
solver$execute_step(T_stop, early_stop = TRUE)
```

so, with `early_stop = TRUE`, the solver terminates early once the number of dummy variables in the active set reaches $T_{\text{stop}}$.

---

## Demo structure

### 1. Basic T-OLS with early stopping

The first part runs T-OLS with **internal normalization** and an intercept on synthetic data, in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part applies external L2 normalization via `LpNormScaler` to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-OLS with `normalize = FALSE` and `intercept = FALSE`.

This demonstrates the externally normalized workflow and allows it to be compared conceptually with the internally normalized workflow from the first part.

### 3. Serialization and warm start

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to `tools_checkpoint.bin` in the demo folder, reloads it, and continues from the checkpoint to the same final stopping level.

Unlike the cpp version, which uses a static `load(filename, X, D)` factory, the R binding hydrates an existing solver instance: a fresh solver is constructed with the same `X`, `D`, and `y`, and the checkpoint is loaded via `solver$load(path)`.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, and equality of the recorded action path. The checkpoint file is removed afterward.

### 4. Memory-mapped data

The fourth part writes the design matrix $X$ to disk via `convert_to_memory_mapped()` and runs T-OLS on the resulting `mmap_matrix` view. The R binding accepts a memory-mapped $X$ while $D$ and $y$ stay in memory; the cpp demo maps all three.

The cpp counterpart generates $p = 500{,}000$ predictors directly on disk; the R port scales this down to $p = 5{,}000$ (with $n = 1{,}000$ and $10p$ dummies). The file `demo_tools_X.bin` is written into the demo folder and removed at the end.

---

## Running

Run from the repository root:

```bash
Rscript R/tsolvers/demo_ts_11_tools/demo_ts_11_tools.R
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

Compared with T-OMP, T-OLS uses a different selection rule: it chooses the candidate that yields the greatest least-squares reduction in residual sum of squares, rather than the one with the largest raw residual correlation. This can make individual steps more expensive, but it may also produce different and sometimes more deliberate active-set decisions.

---

## Notes

- Synthetic data come from the suite-level helper `ts_demo_utils.R::gen_synthetic_data()`, an R replication of the cpp `datagen::SyntheticData` generator (SNR-calibrated Gaussian noise, seed 42), so the selected path can be compared against ground truth.
- Selection diagnostics come from `print_selection` / `print_selection_quality`, which compute TPP/FDP/precision/recall/F1 via the TRexSelectorNeo `compute_*` helpers. The active-dummy count is reconstructed from the signed action path because the R binding's `get_actives()` returns predictors only.
- The serialization example writes a temporary checkpoint file named `tools_checkpoint.bin` and removes it afterward.
- The memory-mapped example writes a temporary file for $X$ only, runs T-OLS on the mapped view, and removes the file at the end.

---

**Last updated**: 2026-07-08
