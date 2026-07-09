# Demo 01: Single-Run T-Rex Selector

## Purpose

Demonstrate basic T-Rex selector usage and output format via two contrasting scenarios:
1. **High-dimensional** ($p > n$): Moderate sample, many predictors
2. **Low-dimensional** ($n > p$): Large sample, moderate predictors

This demo verifies correct API usage and produces human-readable console/file output for inspection. `main()` runs the high-dimensional scenario first, then the low-dimensional one; both use fixed $+1$ coefficients (`rnd_coef = false`).

---

## Data Generation Parameters

Both scenarios share the same true support and coefficients; only the dimensions differ.

### Scenario A: High-Dimensional ($p > n$)
- $n = 150$ observations
- $p = 300$ features
- True support: $\mathcal{S}^* = \{27, 149, 43, 128, 42, 4\}$ ($s = 6$ true actives)
- True coefficients: fixed $\beta_j = 1$ for all active variables
- SNR: $1.0$
- Model: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, \sigma^2 I_n)$
- Data seed: `58`

### Scenario B: Low-Dimensional ($n > p$)
- $n = 5000$ observations
- $p = 1000$ features
- True support: Same as scenario A
- True coefficients: Same as scenario A
- SNR: $1.0$

**Coefficient note**: the `rnd_coef` code path does **not** draw from a uniform distribution — it selects a fixed hardcoded vector $\{-0.4, -0.25, -0.8, 1.1, 2.5, -1.2\}$. Both runs in `main()` use `rnd_coef = false`, i.e. all coefficients equal $1$.

---

## Control Parameters

```
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Max L = 10 × p
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
dummy_distribution = Normal      # Dummy predictors drawn from N(0,1)
lloop_strategy = HCONCAT         # Horizontally concatenated dummy columns
tloop_stagnation_stop = false    # No stagnation early-exit
tloop_max_stagnant_steps = 5     # (inactive: stagnation stop disabled)
parallel_rnd_experiments = false # Sequential dummy experiments
solver_type = TLARS              # Base solver
tFDR = 0.1                       # Target FDR control level
```

---

## Output Files

Two output files are generated in `simulation_results/`, one per scenario. The file stem is built as `d01_trex_basic_n{n}_p{p}.txt`:

### `d01_trex_basic_n150_p300.txt`
**Scenario A output** (high-dimensional case)

The program prints (both to console and file), in order: the selected indices, FDP, TPP, and the internal T-Rex matrices. Actual contents:
```
High-dimensional (p > n)
Generating synthetic data...
Creating maps of data...
Creating T-Rex Selector instance...
Executing T-Rex Selector...
Selected indices: 27 43 128
False Discovery Proportion (FDP): 0.0000
True Positive Proportion (TPP):   0.5000

Adjusted Relative Occurrences (Phi_prime):
...

Phi Matrix (Phi):
...

Estimated FDP Matrix (FDP_hat):
...

R Matrix (R):
...

Voting Grid:
 0.5 0.55  0.6 0.65  0.7 0.75  0.8 0.85  0.9 0.95
```

### `d01_trex_basic_n5000_p1000.txt`
**Scenario B output** (low-dimensional case), same block structure.

There is no timing block and no "number of selected variables / false positives" summary — the program reports the selected index list plus FDP/TPP and the four diagnostic matrices (Phi', Phi, FDP_hat, R) and the voting grid.

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex/demo_trex_01_single_run/demo_trex_01_single_run
```

Output appears on the console and is simultaneously written to the files above.

---

## Expected Results

- **FDP** is well controlled: the committed high-dimensional run reports $\text{FDP} = 0.0$ (no false selections).
- **TPP** at SNR $= 1.0$ is partial in the harder high-dimensional case (the committed run recovers 3 of the 6 actives, $\text{TPP} = 0.5$); the easier low-dimensional case recovers more.
- The **voting grid** spans $\{0.5, 0.55, \ldots, 0.95\}$ and the **R** / **FDP_hat** matrices show how the relative-occurrence threshold and dummy multiplier $L$ evolve across the T-loop.

Results are from a single random seed and vary with the scenario; the high-dimensional case ($p > n$) is intrinsically harder than the low-dimensional one.

---

## Interpretation

- **Console output** shows real-time progress and the final selection summary.
- **File output** provides a permanent record for the single run.
- This demo is NOT a Monte Carlo study; results are from a single random seed and should not be interpreted as population averages.
- Compare FDR/TPR with Demos 02–04 for empirical performance over many repetitions.

---

**Last updated**: 2026-07-08
