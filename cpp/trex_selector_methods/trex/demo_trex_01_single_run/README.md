# Demo 01: Single-Run T-Rex Selector

## Purpose

Demonstrate basic T-Rex selector usage and output format via two contrasting scenarios:
1. **Low-dimensional** ($n > p$): Large sample, moderate predictors
2. **High-dimensional** ($p > n$): Moderate sample, many predictors

This demo verifies correct API usage and produces human-readable console/file output for inspection.

---

## Data Generation Parameters

### Scenario A: Low-Dimensional ($n > p$)
- $n = 5000$ observations
- $p = 1000$ features
- True support: $\mathcal{S}^* = \{27, 149, 43, 128, 42, 4\}$ ($s = 6$ true actives)
- True coefficients: Uniform $\pm 1$ or random $\sim U(-2.5, 2.5)$
- SNR: $1.0$
- Model: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, I_n)$

### Scenario B: High-Dimensional ($p > n$)
- $n = 150$ observations
- $p = 300$ features
- True support: Same as scenario A
- True coefficients: Same as scenario A
- SNR: $1.0$

---

## Control Parameters

```
K = 20                           # Dummy matrix experiments per T-loop iteration
max_dummy_multiplier = 10        # Max L = 10 × p
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
tloop_stagnation_stop = false    # No early exit (single run, let full convergence)
opt_threshold = 0.75             # Optimization grid point
```

---

## Output Files

Two output files are generated in `simulation_results/`:

### `d01_trex_basic_n150_p300.txt`
**Scenario B output** (high-dimensional case)

Example contents:
```
High-dimensional (p > n)
Generating synthetic data...
Creating maps of data...
Setup Control Structures
...
Selection Results
Solver: TLARS
Number of selected variables: 8
Selected variables (0-indexed): [27, 149, 43, 128, 42, 4, ...]
True positives: 6 / 6
False positives: 2
FDR estimate: 0.25
```

### `d01_trex_basic_n5000_p1000.txt`
**Scenario A output** (low-dimensional case)

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_01_single_run
```

Output appears on console and is simultaneously written to the files above.

---

## Expected Results

Both scenarios should produce:
- **Selection size** around 7–10 variables (controls near tFDR = 0.1)
- **FDR** well-controlled (near 0.1, capped by tFDR)
- **True positive count** = 6 (recovery of all true actives) at SNR = 1.0
- **False positives** 1–3 on average (FDR × selected count)

The high-dimensional case ($p > n$) typically requires more dummies ($L$ larger) than low-dimensional due to the increased noise space.

---

## Interpretation

- **Console output** shows real-time progress and final selection summary
- **File output** provides a permanent record for the single run
- This demo is NOT a Monte Carlo study; results are from a single random seed and should not be interpreted as population averages
- Compare FDR/TPR with Demos 02–04 for empirical performance over many repetitions

---

**Last updated**: 2026-07-01
