# Demo 05: Dummy Distribution Comparison

## Purpose

Compare the **dummy distributions** of the T-Rex Selector — the distribution from which the injected dummy (null) variables are drawn — across **three base solver families** (TLARS, TOMP, TAFS), with the L-loop strategy fixed at **STANDARD**. The T-Rex FDR calibration rests on exchangeability between real null predictors and the dummies; since dummies are centered and L2-normalized before entering the solver, their *scale* is immaterial. What this demo probes is whether the *shape* of the dummy distribution (tail weight, discreteness, skewness, sparsity, cross-column dependence) affects the achieved FDR / TPR and the calibrated dummy multiplier $L$ / stopping time $T$ — and whether any such effect depends on the solver family (equiangular LARS path vs greedy forward selection).

---

## Data Generation Parameters

- **Sample size**: $n = 300$
- **Number of features**: $p = 1000$
- **True support**: a random subset of size $s = 10$, **redrawn per trial** (`block_support = false`) by shuffling $\{0, \ldots, p-1\}$ with `std::mt19937(seed + 500000)`, keeping 10 indices and sorting them. A contiguous block support $\{0, \ldots, s-1\}$ is available via the `block_support = true` call, which is compiled but guarded by `if (false)`.
- **True coefficients**: fixed $\beta_j = 1$ (`rnd_coef = false`)
- **SNR grid**: 21 values — $\{0.1, 0.2, \ldots, 2.0\}$ plus $5.0$
- **Monte Carlo repetitions**: `num_MC = 10` trials per solver × distribution × SNR level
- **DGP**: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, \sigma^2 I_n)$, Normal predictors and Normal noise
- **Base solvers** (outer sweep, one result file pair each):
  - **TLARS** — equiangular LARS path (T-loop stagnation stop AUTO-resolves to *disabled*),
  - **TOMP** — greedy orthogonal matching pursuit (stagnation stop AUTO-resolves to *enabled*),
  - **TAFS** — greedy adaptive forward selection with `rho_afs = 1.0` (stagnation stop AUTO-resolves to *enabled*)
- **L-loop strategy**: STANDARD (fixed; fresh i.i.d. dummy matrices at each L-loop iteration)
- **tFDR**: $0.1$

---

## Dummy Distributions Tested

`make_dummy_distributions()` sweeps the distribution types defined by the library `dummygen::Distribution` factory (see `utils/datageneration/utils_dummygen.hpp`), with Student's t evaluated at two tail weights — **twelve rows** in total:

1. **Normal** — $N(0, 1)$; the reference choice.
2. **Uniform** — $U(-\sqrt{3}, \sqrt{3})$; bounded, light tails.
3. **Rademacher** — $\{-1, +1\}$ equiprobable; discrete two-point.
4. **StudentT_df3** — heavy tails (df = 3, unit-variance scaled).
5. **StudentT_df5** — moderately heavy tails (df = 5, unit-variance scaled).
6. **Laplace** — double exponential; heavier-than-normal tails.
7. **Gumbel** — extreme-value distribution; skewed (centered at 0).
8. **Triangle** — bounded, symmetric; $(-\sqrt{6}, 0, \sqrt{6})$.
9. **Logistic** — tail weight between normal and Laplace.
10. **Mammen** — asymmetric golden-ratio two-point distribution.
11. **SparseRad_s0.1** — constrained sparse Rademacher, sparsity $s = 0.1$ (90% zeros, balanced $\pm 1$ entries).
12. **UnifSphere_d5** — uniform on the unit 5-sphere; consecutive 5-column groups are dependent (norm constraint), a deliberate probe of the exchangeability axis.

**Note on UniformSphere**: the library requests dummies in multiples of $p = 1000$, so the sphere dimension must divide $p$; the factory default `dim = 3` would throw for this configuration, hence `dim = 5`.

All distributions use their unit-variance default parameters where applicable.

---

## Control Parameters

Set per solver × distribution (see the solver / distribution loops):

```
solver_type = <varies>           # TLARS / TOMP / TAFS (outer sweep)
solver_params.rho_afs = 1.0      # TAFS regularization (ignored by TLARS/TOMP)
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Upper bound for the adaptive L-loop
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
dummy_distribution = <varies>    # Swept across the 12 rows (inner sweep)
lloop_strategy = STANDARD        # Fresh i.i.d. dummies per L-loop iteration, fixed
tloop_stagnation_stop = AUTO     # Disabled for TLARS, enabled for TOMP/TAFS
tFDR = 0.1                       # Target FDR control level
```

The MC loop is parallelized with OpenMP (`omp_set_num_threads(6)`).

---

## Output Files

All files are written to `simulation_results/data/`. One `.txt` + `.csv` pair is written **per solver**; the stem encodes the solver and the support scenario (`random_support` for the active run):

### Main Result Files
**`demo_trex_05_dummy_distributions_results_n300_p1000_{tlars,tomp,tafs}_random_support.txt`**

Aligned table (written by the shared `save_and_print_mc_results`) with four metric rows — FDR, TPR, Avg L, Avg T — per distribution across the 21 SNR columns. The row label is the distribution name (Normal, Uniform, …, UnifSphere_d5):

```
======================================================================
=== T-Rex Results (averaged over 10 Monte Carlo runs) ===
======================================================================

Solver         Metric    SNR       0.1       0.2  ...       2.0       5.0
--------------------------------------------------------------------------
Normal         FDR              ...
               TPR              ...
               Avg L            ...
               Avg T            ...

Uniform        ...
...
```

### Tidy-Format CSVs
**`demo_trex_05_dummy_distributions_results_n300_p1000_{tlars,tomp,tafs}_random_support.csv`**

Long/stacked format, header column order **`solver,metric,snr,value`** (here the "solver" column holds the distribution name; the base solver is encoded in the file name), with `FDR`, `TPR`, `AvgL`, `AvgT` rows.

### Figures

`./generate_plots.sh` renders one figure set per solver (overview + grouped static figures and an interactive HTML) into `simulation_results/plots/` via the suite-level [../trex_plt_utils.py](../trex_plt_utils.py).

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex/demo_trex_05_mc_sim_dummy_distributions/demo_trex_05_mc_sim_dummy_distributions
```

---

## Key Questions Addressed

1. **Is the FDR calibration robust to the dummy distribution?**
   - Expected: FDR $\leq$ tFDR for all unit-variance i.i.d. distributions, since dummies are L2-normalized and the calibration only relies on exchangeability with the null predictors.

2. **Do heavy tails, skewness, or discreteness matter after L2 normalization?**
   - Compare StudentT_df3 / Laplace (tails), Gumbel / Mammen (skewness), and Rademacher / SparseRad (discreteness) against the Normal reference.

3. **Does within-group dependence (UnifSphere_d5) degrade the calibration?**
   - The sphere rows deliberately violate the independent-columns assumption; watch for FDR or Avg L deviations relative to the i.i.d. rows.

4. **Is the distribution-robustness solver-independent?**
   - Compare the TLARS (equiangular), TOMP, and TAFS (greedy) file pairs: the dummy calibration argument does not depend on the solver family, so the FDR-equivalence pattern should replicate across all three.

---

## Interpretation Guide

**What to look for:**
- **FDR equivalence**: the i.i.d. rows should all keep FDR $\leq$ tFDR with only MC noise between them.
- **TPR consistency**: minor variation due to randomness and the small MC count (`num_MC = 10`).
- **Avg L / Avg T**: whether the L-loop calibrates a different dummy budget for a given distribution shape.

**Practical significance:**
- If all shapes perform equivalently, the Normal default is justified and the distribution choice is a free parameter (e.g., Rademacher / SparseRad admit cheaper generation or sparse storage).
- Systematic deviations for the dependent UnifSphere rows would illustrate why the exchangeability requirement in `utils_dummygen.hpp` matters.

**Recorded results** (`simulation_results/`, 2026-07-11, num_MC = 10, TLARS + TOMP + TAFS): the FDR-equivalence pattern replicates across all three solver families. For every solver, all twelve distribution rows behave equivalently — the per-distribution mean FDR spread is tiny (TLARS 0.046–0.057, TOMP 0.025–0.035, TAFS 0.024–0.038), TPR curves are statistically indistinguishable and reach 1.0 at SNR 5.0, and even the dependent UnifSphere_d5 rows show no systematic deviation (they yield the *lowest* mean FDR under TOMP) — at this ($n$, $p$, dim) the within-group sphere dependence is too weak to break the calibration after L2 normalization. Isolated FDR exceedances up to ≈0.2 at SNR = 0.2 are 10-trial MC noise on near-empty selections and appear for all rows alike (mirroring Demo 04); from SNR ≥ 0.3 the greedy solvers stay ≤ 0.07 and TLARS ≤ 0.12. The greedy TOMP/TAFS runs (stagnation stop enabled) are systematically more conservative than TLARS (mean FDR ≈ 0.03 vs ≈ 0.05) with a tighter stopping-time range (Avg T ∈ [4, 9] vs [2, 23]).

---

**Last updated**: 2026-07-11
