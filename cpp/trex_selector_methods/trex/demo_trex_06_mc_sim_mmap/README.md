# Demo 06: Monte Carlo Simulation with Memory-Mapped Data

## Purpose

Validate **T-Rex selector performance under memory-mapped (mmap) I/O** over many Monte Carlo trials, mirroring the single-run patterns of Demo 05. Two scenarios are run: (C) in-memory $\mathbf{X}$ with mmap dummy matrices $\mathbf{D}$ + solver serialization, and (D) a fully memory-mapped pipeline where each MC iteration creates and destroys its own mmap-backed $\mathbf{X}/\mathbf{y}$ files.

---

## Data Generation Parameters

Both scenarios use the same high-dimensional configuration:

- **Sample size**: $n = 300$
- **Number of features**: $p = 1000$
- **True support**: a set of **10 unique random indices** drawn once from `std::mt19937 rng(24)` (uniform over $\{0, \ldots, p-1\}$), fixed across all MC trials
- **True coefficients**: fixed $\beta_j = 1$ (`rnd_coef = false`)
- **SNR range**: $\{0.1, 0.5, 1.0, 2.0, 5.0\}$
- **Monte Carlo repetitions**:
  - **Demo C**: `num_MC = 10`
  - **Demo D**: `num_MC = 500`
- **Base solver**: TLARS only
- **DGP**: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, \sigma^2 I_n)$
- **tFDR**: $0.1$

---

## Scenarios

### Demo C: MC — In-Memory $\mathbf{X}$ + `use_memory_mapping = true`
- $\mathbf{X}$ generated in RAM per trial; dummy matrices $\mathbf{D}$ memory-mapped and solver LARS-path state serialized to disk.
- Purpose: verify the D-mmap + solver-serialization pipeline yields stable, reproducible FDR/TPR over many runs.

### Demo D: MC — Fully Memory-Mapped Pipeline ($\mathbf{X}$ + $\mathbf{D}$ + solver serialization)
- Each MC iteration writes its own mmap-backed $\mathbf{X}/\mathbf{y}$ files (per-thread paths derived from the filestem `X_mmap_mc` and `omp_get_thread_num()`), which are removed by an RAII `MmapFileGuard` at the end of the iteration scope — exception-safe.

---

## Control Parameters

```
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Max L = 10p
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
dummy_distribution = Normal      # Dummy predictors drawn from N(0,1)
lloop_strategy = HCONCAT         # Horizontally concatenated dummy columns
tloop_stagnation_stop = true     # Early exit when R_mat stagnates
tloop_max_stagnant_steps = 7     # Stagnation window (7) → "sw7" in the output stem
opt_threshold = 0.75             # Optimization grid point
parallel_rnd_experiments = false # Sequential dummy experiments
use_memory_mapping = true        # D mmap + solver serialization
solver_type = TLARS              # Base solver
tFDR = 0.1                       # Target FDR control level
```

The MC loop is parallelized with OpenMP (`omp_set_num_threads(6)`).

---

## Output Files

Written to `simulation_results/`. Each scenario writes one `.txt` and one `.csv`, with stems encoding `n`, `p`, and the stagnation window (`sw7`):

- **Demo C**: `d03_trex_mmap_demo_c_n300_p1000_sw7.{txt,csv}`
- **Demo D**: `d03_trex_mmap_demo_d_n300_p1000_sw7.{txt,csv}`

The `.txt` file (written by the shared `save_and_print_mc_results`) is an aligned table with FDR, TPR, Avg L, and Avg T rows for the single TLARS solver across the 5 SNR columns:

```
======================================================================
=== T-Rex Results (averaged over 10 Monte Carlo runs) ===
======================================================================

Solver         Metric    SNR       0.1       0.5       1.0       2.0       5.0
------------------------------------------------------------------------------
TLARS          FDR              ...
               TPR              ...
               Avg L            ...
               Avg T            ...
```

(Demo D's header reports "averaged over 500 Monte Carlo runs".) The `.csv` uses the tidy header **`solver,metric,snr,value`** with `FDR`, `TPR`, `AvgL`, `AvgT` rows.

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex/demo_trex_06_mc_sim_mmap/demo_trex_06_mc_sim_mmap
```

`main()` runs Demo C (`num_MC = 10`) then Demo D (`num_MC = 500`).

---

## Key Questions Addressed

1. **Does FDR control hold with mmap storage?**
   - Expected: FDR $\leq$ tFDR — the storage medium is transparent to the algorithm.

2. **Does the fully-mmapped pipeline (Demo D) manage its per-iteration files safely?**
   - Expected: yes; RAII `MmapFileGuard` cleanup runs even if `select()` throws.

3. **Are results stable across many trials?**
   - Demo D's 500 trials give tighter FDR/TPR estimates than Demo C's 10.

---

## Interpretation Guide

**What to look for:**
- **FDR stability**: FDR should stay $\leq$ tFDR across SNR.
- **TPR progression**: increases with SNR toward full recovery.
- **Avg L / Avg T**: dummy-multiplier and stopping-time behaviour under the mmap pipeline.
- **Reliability**: no leaked mmap files, no crashes.

**Comparison with Demo 05 (single-run mmap):**
- Demo 05 validates basic mmap correctness (one run per scenario).
- Demo 06 validates robustness and FDR guarantees over many MC trials.

---

**Last updated**: 2026-07-08
