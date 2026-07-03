# Demo: T-Stepwise Solver

## Purpose

This demo illustrates the **T-Stepwise (Terminating Stepwise Regression)** solver.

T-Stepwise is a greedy variable-selection method that, at each step, adds the variable most correlated with the current residual. The demo shows how the solver is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

---

## Key idea

Like the other terminating solvers, T-Stepwise is run only up to a prescribed stopping level $T_{\text{stop}}$:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TStepwise}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form

```cpp
executeStep(T_stop, /*early_stop=*/true)
```

so the solver terminates once the requested number of steps has been reached.

---

## Demo structure

### 1. Basic T-Stepwise with early stopping

The first part runs T-Stepwise with **internal normalization** and an intercept on synthetic data. It is executed in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part applies external L2 normalization to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-Stepwise with internal normalization disabled.

This demonstrates the externally normalized workflow and allows it to be compared conceptually with the internally normalized workflow used in the first part.

### 3. Serialization and path consistency

The third part runs a reference fit to a final stopping level, then runs a shorter partial fit, saves the solver state to disk, reloads it, and continues from the checkpoint to the same final stopping level.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, and equality of the recorded action path.

### 4. Memory-mapped data

The fourth part generates synthetic data directly on disk, retrieves memory-mapped views of $X$, $D$, and $y$, and runs T-Stepwise on those mapped objects.

This shows that the solver can also be used in large-scale workflows where the matrices are stored on disk rather than held entirely in ordinary in-memory arrays.

---

## Running

```bash
./build/debug/bin/tsolvers/demo_ts_03_tstepwise/demo_ts_03_tstepwise
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

Because T-Stepwise uses a greedy residual-correlation rule rather than an equiangular LARS-style update, its selection order may differ noticeably from T-LARS on the same problem.

---

## Notes

- The synthetic data are generated with known support and coefficients, so the selected path can be compared against ground truth.
- The first demo uses a lower signal-to-noise ratio and a slightly different default coefficient pattern than the external-normalization and memory-mapped examples, so behavior across sections is not meant to be numerically identical.
- The serialization example writes a temporary checkpoint file named `tstepwise_checkpoint.bin` and removes it afterward.
- The memory-mapped example writes temporary files for $X$, $D$, and $y$, runs T-Stepwise on those mapped views, and removes those files at the end.

---

**Last updated**: 2026-07-03
