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

A terminating solver does not compute a full path over all possible steps. Instead, it is run only up to a prescribed stopping level $T_{\text{stop}}$:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TLARS}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form:

```cpp
executeStep(T_stop, /*early_stop=*/true)
```

which stop the solver once the requested number of steps has been reached.

---

## Demo structure

### 1. Basic T-LARS with early stopping

The first part runs T-LARS with **internal normalization** and an intercept on synthetic data, and it is executed in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part preprocesses the signal matrix and dummy matrix externally using L2-based column normalization, centers the response manually, and then runs T-LARS with internal normalization disabled.

This section is useful for comparing the solver’s external-normalization workflow against the internally normalized workflow from the first demo.

### 3. Serialization and warm start

The third part runs a reference T-LARS fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues to the same final stopping level.

The reloaded solver is compared against the reference run through the number of steps, the final RSS difference, the final $R^2$ difference, and the equality of the recorded action path.

### 4. Memory-mapped data

The fourth part generates synthetic data directly on disk, retrieves memory-mapped matrix views, and runs T-LARS on those mapped objects.

This demonstrates that the solver can also be used in workflows where the data are too large or too inconvenient to keep entirely in ordinary in-memory matrices.

---

## Running

```bash
./build/debug/bin/tsolvers/demo_ts_01_tlars/demo_ts_01_tlars
```

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

- The demo uses synthetic data with known support and coefficients, so the selected path can be compared against ground truth.
- The in-memory examples use both ordinary dense matrices and externally normalized matrices.
- The serialization example writes a temporary checkpoint file and removes it after the comparison.
- The memory-mapped example writes temporary files for $X$, $D$, and $y$, runs the solver on mapped views, and then removes those files at the end.

---

**Last updated**: 2026-07-03
