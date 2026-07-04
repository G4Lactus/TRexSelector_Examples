# Demo 06: Monte Carlo Simulation with Memory-Mapped Data

## Purpose

Validate **T-Rex selector performance with memory-mapped (mmap) I/O** over multiple Monte Carlo trials. This demo ensures FDR control and power are maintained when working with disk-resident data matrices, not just single-run cases.

---

## Data Generation Parameters

- **Sample size**: $n = 5000$
- **Number of features**: $p = 10000$
- **True support**: $\mathcal{S}^* = \{$first 20 features$\}$ (fixed)
- **True coefficients**: Uniform $\pm 1$
- **SNR range**: $\{0.5, 1.0, 2.0, 5.0\}$
- **Monte Carlo repetitions**: 100 trials per SNR level (fewer than Demos 02–04 due to computational cost)
- **Data storage**: Both $\mathbf{X}$ and $\mathbf{D}$ memory-mapped to disk
- **DGP**: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, I_n)$

---

## Control Parameters

```
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Max L = 10p
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
tloop_stagnation_stop = true     # Early exit
opt_threshold = 0.75             # Optimization grid point
target_fdr = 0.1                 # Target FDR
use_memory_mapping = true        # Fully enable mmap for X and D
parallel_rnd_experiments = false  # Sequential dummy runs
```

---

## Output Files

### Main Result File
**`demo_trex_06_mc_sim_mmap_trex_results_n5000_p10000.txt`**

Tabular format with FDR/TPR by solver and SNR (all solvers with mmap backend):

```
======================================================================
=== T-Rex Results with Memory-Mapped Data (averaged over 100 MC runs) ===
======================================================================

Solver         Metric    SNR       0.5       1.0       2.0       5.0
--------------------------------------------------------------------
TLARS          FDR              0.0412    0.0398    0.0385    0.0271
               TPR              0.2580    0.7340    0.9850    1.0000

TLASSO         FDR              0.0425    0.0415    0.0392    0.0278
               TPR              0.2540    0.7325    0.9840    1.0000

... (additional solvers)
```

### CSV Format
**`demo_trex_06_mc_sim_mmap_trex_results_n5000_p10000.csv`**

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_06_mc_sim_mmap
```

**Computation time**: 1–2 hours (significantly longer than Demos 02–04 due to:
- Larger data ($p = 10000$ vs. $p = 1000$)
- Mmap I/O overhead
- Fewer MC trials (100 vs. 200) to keep runtime manageable

**Disk space**: ~5–10 GB for temporary mmap files (ensure available space in `/tmp/` or configured directory)

---

## Key Questions Addressed

1. **Does FDR control hold with mmap storage?**
   - Expected: Yes, FDR $\leq$ tFDR (storage medium is transparent to algorithm)

2. **Does mmap introduce memory leaks or file handle exhaustion?**
   - Expected: No; RAII cleanup ensures proper resource management

3. **Is performance degradation acceptable for the use case?**
   - Expected: ~1.5–3× slower than in-memory, but necessary for ultra-large datasets

---

## Interpretation Guide

**What to look for:**
- **FDR stability**: Should match Demo 02 results (modulo SNR/support differences)
- **TPR patterns**: Similar progression with SNR as in-memory versions
- **Solver agreement**: All solvers maintain identical FDR/TPR
- **Reliability**: No errors, crashes, or incomplete output files

**Comparison with Demo 05 (single-run mmap):**
- Demo 05 validates basic mmap correctness (one trial per scenario)
- Demo 06 validates robustness over many trials and FDR guarantees

---

## Cleanup

After the demo completes, remove temporary mmap files:
```bash
rm -rf /tmp/trex_mmaps/
```

---

**Last updated**: 2026-07-01
