# Demo: T-LASSO Solver

## Purpose

This demo illustrates the **T-LASSO (Terminating Least Absolute Shrinkage and Selection Operator)** solver.

It shows how T-LASSO is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

In addition to the standard selection diagnostics, the demo also reports T-LASSO-specific path diagnostics such as variable removals and cycling behavior.

---

## Key idea

Like the other terminating solvers, T-LASSO stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set, rather than tracing a full path:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TLASSO}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form

```cpp
executeStep(T_stop, /*early_stop=*/true)
```

where `T_stop` is the number of dummy variables that must enter the active set before the solver terminates early: with `early_stop=true` the run stops as soon as the count of active dummies reaches $T_{\text{stop}}$.

---

## Demo structure

### 1. Basic T-LASSO with early stopping

The first part runs T-LASSO with **internal normalization** and an intercept on synthetic data. It is executed in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables, selection-quality diagnostics, and T-LASSO-specific quantities such as the number of removals and the cycling ratio.

### 2. External normalization

The second part applies external L2 normalization to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-LASSO with internal normalization disabled.

This section demonstrates the externally normalized workflow and again reports the resulting selection path, selection quality, removals, and cycling statistics.

### 3. Serialization and path consistency

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues from the checkpoint to the same final stopping level.

The reloaded path is compared against the reference run using the final number of steps, the final RSS difference, the final $R^2$ difference, equality of the recorded action path, and the T-LASSO-specific removal and cycling diagnostics.

### 4. Memory-mapped data

The fourth part generates synthetic data directly on disk, retrieves memory-mapped views of $X$, $D$, and $y$, and runs T-LASSO on those mapped objects.

This shows that the solver can also be used in workflows where large matrices are stored on disk rather than held entirely in standard in-memory matrices.

---

## Running

```bash
./build/debug/bin/tsolvers/demo_ts_02_tlasso/demo_ts_02_tlasso
```

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

- The synthetic data are generated with known support and coefficients, so the selection path can be compared against ground truth.
- The demo uses both in-memory and memory-mapped workflows.
- The serialization example writes a temporary checkpoint file named `tlasso_checkpoint.bin` and removes it afterward.
- The memory-mapped example writes temporary files for $X$, $D$, and $y$, runs T-LASSO on those mapped views, and removes those files at the end.
- For a dedicated cross-check against related LARS/LASSO reference behavior, see the validation folder described in `validation_ts_02_tlars_tlasso_rcompare`.

---

**Last updated**: 2026-07-08
