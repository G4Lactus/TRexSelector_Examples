# Demo 03: Monte Carlo Simulation with Variable Support Size

## Purpose

Investigate how T-Rex selector performance (FDR and TPR) varies as **sparsity level** changes. This demo shows robustness across different support sizes and tests whether FDR control is maintained when the true signal strength ($s$) varies.

---

## Data Generation Parameters

- **Sample size**: $n = 300$
- **Number of features**: $p = 1000$
- **True support size**: Variable across Monte Carlo trials
  - Trials may have different random subsets of varying cardinalities
  - Typical range: $s \in \{5, 10, 15, 20, \ldots\}$ (configurable)
- **True coefficients**: Uniform $\pm 1$ or random $\sim U(-2.5, 2.5)$
- **SNR grid**: Multiple levels (e.g., $\{0.5, 1.0, 2.0, 5.0\}$)
- **Monte Carlo repetitions**: 200 trials per SNR level
- **DGP**: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, I_n)$

---

## Control Parameters

```
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Max dummies L = 10p
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
tloop_stagnation_stop = true     # Early exit when R_mat stagnates
tloop_max_stagnant_steps = 7     # Stagnation window
opt_threshold = 0.75             # Optimization grid point
target_fdr = 0.1                 # Target FDR
```

---

## Output Files

Similar structure to Demo 02:

### Main Result File
**`demo_trex_03_mc_sim_variable_support_trex_results_*.txt`**

Tabular FDR and TPR by solver and SNR (when support is variable, results are aggregated across trials).

### CSV Format
**`demo_trex_03_mc_sim_variable_support_trex_results_*.csv`**

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_03_mc_sim_variable_support
```

**Computation time**: Similar to Demo 02 ($\approx 10–30$ minutes).

---

## Key Questions Addressed

1. **Does FDR control degrade when support size varies?**
   - Expected: FDR remains controlled (algorithm is support-size-agnostic)

2. **How does TPR change with sparsity?**
   - Expected: More sparse signals (smaller $s$) may have higher relative power per variable

3. **Is there interaction between SNR and support size?**
   - Expected: Weak signals with small $s$ are harder to detect; strong signals are easier regardless of $s$

---

## Interpretation Guide

**What to look for:**
- **FDR stability**: Should remain $\leq$ tFDR even for variable support sizes
- **TPR patterns**: Compare against fixed-support Demo 02; larger variations in $s$ may introduce TPR variability
- **Solver consistency**: All solvers should maintain agreement regardless of support size

**Comparison with Demo 02:**
- Demo 02 (fixed $s = 10$) provides a control baseline
- Demo 03 (variable $s$) shows robustness under realistic uncertainty about true sparsity

---

**Last updated**: 2026-07-01
