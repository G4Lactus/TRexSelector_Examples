# Demo: T-AFS (Adaptive Forward Stepwise) Solver

## Purpose

This demo illustrates the **T-AFS (Terminating Adaptive Forward Stepwise)** solver.

T-AFS is a forward-selection method with an adaptive decision rule that adjusts the variable-selection criterion as the active set evolves. The demo shows how the solver is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

---

## Key idea

Like the other terminating solvers, T-AFS is run only up to a prescribed stopping level $T_{\text{stop}}$:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TAFS}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form

```cpp
executeStep(T_stop, /*early_stop=*/true)
```

so the solver terminates once the requested number of steps has been reached.

---

## Demo structure

### 1. Basic T-AFS with early stopping

The first part runs T-AFS with **internal normalization** and an intercept on synthetic data, in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part applies external L2 normalization to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-AFS with internal normalization disabled.

This demonstrates the externally normalized workflow and allows it to be compared conceptually with the internally normalized workflow from the first part.

### 3. Serialization and warm start

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues from the checkpoint to the same final stopping level.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, and equality of the recorded action path.

### 4. Memory-mapped data

The fourth part generates synthetic data directly on disk, retrieves memory-mapped views of $X$, $D$, and $y$, and runs T-AFS on those mapped objects.

This shows that the solver can also be used in workflows where large matrices are stored on disk rather than held entirely in ordinary in-memory matrices.

---

## Running

```bash
./build/debug/bin/tsolvers/demo_ts_12_tafs/demo_ts_12_tafs
```

The demo prints its diagnostics to the console, including configuration summaries, selected variables, selection-quality measures, serialization checks, and memory-mapped execution status messages.

---

## What to look for

When reading the console output, focus on the following points:

- whether the selected variables recover the known true support reasonably well,
- how the solver behaves under both high-dimensional and low-dimensional settings,
- whether externally normalized and internally normalized workflows both execute successfully,
- whether the serialized solver state can be reloaded and continued without changing the final path,
- whether the memory-mapped workflow completes successfully and cleans up its temporary files.

Conceptually, T-AFS is a useful counterpart to plain T-Stepwise because its adaptive selection rule may behave differently once the active set begins to grow, especially in correlated designs.

---

## Notes

- The synthetic data are generated with known support and coefficients, so the selected path can be compared against ground truth.
- The serialization example writes a temporary checkpoint file and removes it afterward.
- The memory-mapped example writes temporary files for $X$, $D$, and $y$, runs T-AFS on those mapped views, and removes those files at the end.

---

**Last updated**: 2026-07-03
