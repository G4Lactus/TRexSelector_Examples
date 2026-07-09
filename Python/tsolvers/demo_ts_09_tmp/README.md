# Demo: T-MP (Matching Pursuit) Solver (Python)

## Purpose

This demo illustrates the **T-MP (Terminating Matching Pursuit)** solver.

T-MP is the simplest pursuit-style greedy method in this group: it repeatedly selects a variable based on residual correlation and updates coefficients incrementally without performing a full orthogonal refit on the active set. The demo shows how the solver is used in four settings:

1. basic execution with early stopping,
2. external versus internal normalization,
3. serialization and warm-start continuation,
4. operation on memory-mapped data.

It mirrors `cpp/tsolvers/demo_ts_09_tmp`.

---

## Key idea

Like the other terminating solvers, T-MP stops early once a prescribed number $T_{\text{stop}}$ of dummy variables has entered the active set:

$$
\hat{\boldsymbol{\beta}}^{(T)} = \mathrm{TMP}(\mathbf{X}, \mathbf{D}, \mathbf{y}; T_{\text{stop}}).
$$

In the demo, this is done through calls of the form

```python
solver.executeStep(T_stop, early_stop=True)
```

so, with `early_stop=True`, the solver terminates early once the number of dummy variables in the active set reaches $T_{\text{stop}}$.

The solver class is `TMP_Solver` from `trex_selector_neo.tsolvers.omp_based`.

---

## Demo structure

### 1. Basic T-MP with early stopping

The first part runs T-MP with **internal normalization** and an intercept on synthetic data, in both a high-dimensional regime $(p > n)$ and a low-dimensional regime $(n > p)$.

After execution, the demo prints the selected variables and selection-quality diagnostics relative to the known true support.

### 2. External normalization

The second part applies external L2 normalization (`LpNormScaler` with `NormType.L2`) to both the signal matrix and the dummy matrix, centers the response manually, and then runs T-MP with internal normalization disabled.

This demonstrates the externally normalized workflow and allows it to be compared conceptually with the internally normalized workflow from the first part.

### 3. Serialization and warm start

The third part runs a reference fit to a final stopping level, then separately runs a shorter partial fit, saves the solver state to disk, reloads it, and continues from the checkpoint to the same final stopping level.

Unlike cpp, which uses a static `load(filename, X, D)` factory, the Python binding hydrates an existing solver instance: a new `TMP_Solver` is constructed with the same `X`, `D`, and `y`, and then `solver.load(path)` restores the checkpointed state.

The reloaded run is compared against the reference path using the final number of steps, the final RSS difference, the final $R^2$ difference, and equality of the recorded action path.

### 4. Memory-mapped data

The fourth part writes $X$, $D$, and $y$ to disk with `numpy_to_memmap()` and runs T-MP on the zero-copy mmap views (via `.to_numpy()`).

Note: the cpp demo generates $p = 500{,}000$ on disk; the Python port scales this down to $p = 5{,}000$ (with $n = 1000$ and $10p$ dummies).

---

## Running

The demo is a script — run it from the repo root with:

```bash
.venv/bin/python Python/tsolvers/demo_ts_09_tmp/demo_ts_09_tmp.py
```

There is no build step; it requires the `trex_selector_neo` package installed in the repo `.venv`.

The demo prints its diagnostics to the console, including configuration summaries, selected variables, selection-quality measures, serialization checks, and memory-mapped execution status messages.

---

## What to look for

When reading the console output, focus on the following points:

- whether the selected variables recover the known true support reasonably well,
- how the solver behaves under both high-dimensional and low-dimensional settings,
- whether externally normalized and internally normalized workflows both execute successfully,
- whether the serialized solver state can be reloaded and continued without changing the final path,
- whether the memory-mapped workflow completes successfully and cleans up its temporary files.

Compared with T-OMP, T-MP does not perform a full orthogonal refit on the active set, so it can converge more slowly and may revisit the same variable multiple times. It therefore serves as a useful baseline for comparison against more elaborate pursuit variants.

---

## Notes

- Synthetic data come from the suite-level helper `ts_demo_utils.gen_synthetic_data()` (numpy replication of the cpp `datagen::SyntheticData`, SNR-calibrated Gaussian noise, seed 42), so the selected path can be compared against ground truth; selection diagnostics use `print_selection` / `print_selection_quality` (TPP/FDP/precision/recall/F1 via `trex_selector_neo.utils`), with 0-based support indices matching cpp.
- The thread count is set to 6 via `trex_selector_neo.utils.set_num_threads`, mirroring the cpp OpenMP setup.
- The serialization example writes a temporary checkpoint file `tmp_checkpoint.bin` into the demo folder and removes it afterward.
- The memory-mapped example writes `demo_tmp_{X,D,y}.bin` into the demo folder, runs T-MP on those mapped views, and removes those files at the end.

---

**Last updated**: 2026-07-08
