# Demo 04: L-Loop Calibration Strategies Comparison

## Purpose

Compare different **L-loop calibration strategies** to understand how the method for selecting the number of dummy variables $L$ affects FDR control and computational efficiency. This demo validates that the algorithm is robust to the L-calibration choice.

---

## Data Generation Parameters

- **Sample size**: $n = 300$
- **Number of features**: $p = 1000$
- **True support**: Random subset of size $s = 10$ (varies across MC trials)
- **True coefficients**: Uniform $\pm 1$
- **SNR range**: $\{0.1, 1.0, 5.0\}$
- **Monte Carlo repetitions**: 200 trials per SNR level
- **DGP**: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, I_n)$

---

## L-Loop Strategies Tested

The demo compares different approaches to the L-loop, which calibrates the number of dummy variables:

1. **Standard (STANDARD)**: Increment $L$ sequentially ($L = 1, 2, 3, \ldots$) until FDP estimate stabilizes
2. **Doubling**: Exponential growth ($L = 1, 2, 4, 8, \ldots$) for faster convergence
3. **Fixed (FIXED)**: Use a pre-computed $L$ value (e.g., $L = \max(2, \lfloor p / 20 \rfloor)$)
4. Other variants (algorithm-specific)

---

## Control Parameters

```
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Max L = 10p
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
tloop_stagnation_stop = true     # Early exit when R_mat stagnates
tloop_max_stagnant_steps = 7     # Stagnation window
opt_threshold = 0.75             # Optimization grid point
target_fdr = 0.1                 # Target FDR
lloop_strategy = VARIES          # Main variable: switch between strategies
```

---

## Output Files

### Main Result File
**`demo_trex_04_mc_sim_lloop_strategies_trex_results_*.txt`**

Format shows FDR/TPR columns for different L-loop strategies (instead of different solvers):

```
======================================================================
=== T-Rex Results by L-Loop Strategy (averaged over 200 MC runs) ===
======================================================================

Strategy       Metric    SNR       0.1       1.0       5.0
-----------------------------------------------------------
STANDARD       FDR              0.0050    0.0427    0.0279
               TPR              0.0030    0.7360    1.0000

DOUBLING       FDR              0.0055    0.0431    0.0282
               TPR              0.0025    0.7340    0.9995

FIXED          FDR              0.0045    0.0425    0.0275
               TPR              0.0028    0.7355    1.0000

...
```

### CSV Format
**`demo_trex_04_mc_sim_lloop_strategies_trex_results_*.csv`**

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_04_mc_sim_lloop_strategies
```

**Computation time**: Approximately 10–30 minutes.

---

## Key Questions Addressed

1. **Do all L-loop strategies maintain FDR control?**
   - Expected: Yes, despite different calibration paths

2. **Is there a trade-off between FDR control and computational speed?**
   - Expected: Some strategies converge faster (fewer L evaluations) but may require more total dummy runs

3. **Which strategy is most efficient?**
   - Expected: Exponential (doubling) or fixed strategies may be faster than STANDARD

---

## Interpretation Guide

**What to look for:**
- **FDR equivalence**: All strategies should produce FDR $\leq$ tFDR
- **TPR consistency**: Minor variations due to randomness; should be statistically indistinguishable
- **Computational overhead**: Compare wall-clock time across strategies (visible in console output)

**Practical significance:**
- If all strategies are equivalent in FDR/TPR, the fastest one is preferred
- L-loop choice impacts runtime but not selection quality (under proper control)

---

**Last updated**: 2026-07-01
