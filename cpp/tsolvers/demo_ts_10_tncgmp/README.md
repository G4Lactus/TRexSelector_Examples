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

Like the other terminating solvers, T-NCGMP is run only up to a prescribed stopping level $T_{\text{stop}}$:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TNCGMP}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form

```cpp
executeStep(T_stop, /*early_stop=*/true)
```

so the solver terminates once the requested number of steps has been reached.

---

## Demo structure

### 1. Basic T-NCGMP with early stopping

The first part runs T-NCGMP with **internal normalization** and an intercept on synthetic data, in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part applies external L2 normalization to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-NCGMP with the same solver variant on the transformed data.

This demonstrates the externally normalized workflow and allows it to be compared conceptually with the internally normalized workflow from the first part.

### 3. Serialization and warm start

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues from the checkpoint to the same final stopping level.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, and equality of the recorded action path.

### 4. Memory-mapped data

The fourth part generates synthetic data directly on disk, retrieves memory-mapped views of $X$, $D$, and $y$, and runs T-NCGMP on those mapped objects.

This shows that the solver can also be used in workflows where large matrices are stored on disk rather than held entirely in ordinary in-memory matrices.

---

## Running

```bash
./build/debug/bin/tsolvers/demo_ts_10_tncgmp/demo_ts_10_tncgmp
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

Conceptually, T-NCGMP extends the simpler pursuit variants with a norm-corrective step, so it is a natural place to watch for improved stability and more controlled long-path behavior relative to plain T-MP.

---

## Notes

- The synthetic data are generated with known support and coefficients, so the selected path can be compared against ground truth.
- The serialization example writes a temporary checkpoint file named `tncgmp_checkpoint.bin` and removes it afterward.
- The memory-mapped example writes temporary files for $X$, $D$, and $y$, runs T-NCGMP on those mapped views, and removes those files at the end.

---

**Last updated**: 2026-07-03
