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

Like the other terminating solvers, T-OLS is run only up to a prescribed stopping level $T_{\text{stop}}$:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TOOLS}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form

```cpp
executeStep(T_stop, /*early_stop=*/true)
```

so the solver terminates once the requested number of steps has been reached.

---

## Demo structure

### 1. Basic T-OLS with early stopping

The first part runs T-OLS with **internal normalization** and an intercept on synthetic data, in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part applies external L2 normalization to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-OLS with internal normalization disabled.

This demonstrates the externally normalized workflow and allows it to be compared conceptually with the internally normalized workflow from the first part.

### 3. Serialization and warm start

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues from the checkpoint to the same final stopping level.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, and equality of the recorded action path.

### 4. Memory-mapped data

The fourth part generates synthetic data directly on disk, retrieves memory-mapped views of $X$, $D$, and $y$, and runs T-OLS on those mapped objects.

This shows that the solver can also be used in workflows where large matrices are stored on disk rather than held entirely in ordinary in-memory matrices.

---

## Running

```bash
./build/debug/bin/tsolvers/demo_ts_11_tools/demo_ts_11_tools
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

Compared with T-OMP, T-OLS uses a different selection rule: it chooses the candidate that yields the greatest least-squares reduction in residual sum of squares, rather than the one with the largest raw residual correlation. This can make individual steps more expensive, but it may also produce different and sometimes more deliberate active-set decisions.

---

## Notes

- The synthetic data are generated with known support and coefficients, so the selected path can be compared against ground truth.
- The serialization example writes a temporary checkpoint file named `tools_checkpoint.bin` and removes it afterward.
- The memory-mapped example writes temporary files for $X$, $D$, and $y$, runs T-OLS on those mapped views, and removes those files at the end.

---

**Last updated**: 2026-07-03
