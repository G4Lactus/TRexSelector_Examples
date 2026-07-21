# Demo 03: T-Rex+GVS on Mixed Blocks (Python)

## Purpose

Evaluate T-Rex+GVS on a more realistic grouped design: four contiguous equicorrelated blocks of unequal size,
 placed in **random order** with **random noise gaps** between and around them, and including one **inactive
 equicorrelated trap block**.
 The design stresses whether the pipeline can recover the active groups while correctly rejecting a null block that
 shares their correlation structure.
 Selection and evaluation are per variable (see
 [What is actually measured](../README.md#what-is-actually-measured-in-these-demos)).

This is the Python port of the C++ demo
[demo_trex_gvs_03_mc_sim_mixed_blocks](../../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_03_mc_sim_mixed_blocks/),
using `trex_selector_neo`'s `TRexGVSSelector`.

---

## Data Generation Parameters (`dgp_mixed_blocks_snr`)

We consider again the linear model:

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon},
\qquad \boldsymbol{\epsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma_{\varepsilon}^2 \boldsymbol{I}_n)
$$

- $\boldsymbol{y} \in \mathbb{R}^n$ is the response vector.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ is the design matrix.
- $\boldsymbol{\beta} \in \mathbb{R}^p$ is the coefficient vector, with $s$ nonzero entries.
- $\boldsymbol{\epsilon}$ is the noise vector, i.i.d. standard normal.
- $\sigma_{\varepsilon}^2$ is the noise variance, calibrated to achieve a target linear signal-to-noise ratio (SNR).
- $n = 200$, $p = 500$, $s = 150$ (the sum of the three active blocks).

Each block of the design matrix $\boldsymbol{X}$ follows the equicorrelated latent-factor model:

$$
X_{ij} = Z_{i,\,g(j)} + \sigma_x\, \xi_{ij}, \qquad \xi_{ij} \sim \mathcal{N}(0,1),
$$

- Four contiguous blocks of sizes $\{20, 50, 80, 65\}$ (215 grouped columns in total).
- **Blocks 0–2 are active** ($\beta_j = 3$); **block 3 is an inactive trap** ($\beta_j = 0$) carrying the same
   equicorrelated structure as the active blocks.
- Blocks are placed in a **random order** with **random-width noise gaps** (5 gap slots: before, between, and after
   the blocks) filled with i.i.d. $\mathcal{N}(0,1)$ background columns.
- Within-block correlation $\rho = 1/(1 + \sigma_x^2)$, same convention as Demo 01.

---

## Control Parameters

```text
K = 20                          # Random experiments per T-loop iteration
tFDR = 0.1                      # Target FDR
corr_max = 0.98                 # HAC auto-clustering correlation threshold
hc_linkage = Single             # Single-linkage HAC
lambda2_method = CV_1SE_CCD     # Elastic-net penalty selection
MC = 200                        # Monte Carlo repetitions per grid point
```

---

## Methods Compared

Three T-Rex+GVS solver variants: **TENET** (elastic net) [[3]](#references),
**TENET_AUG** (row-augmented elastic net) [[1]](#references), and **TIENET_AUG** (informed
elastic net) [[2]](#references). All use single-linkage HAC grouping
and `CV_1SE_CCD` for the $\lambda_2$ selection.

---

## The Sweep

A single **2-D SNR $\times$ $\rho$ grid**, over
$\mathrm{SNR} \in \{0.2, 0.5, 1, 2, 5\}$ and
$\rho \in \{0.30, 0.50, 0.70, 0.85, 0.90, 0.95, 0.99\}$
(with $\sigma_x = \sqrt{(1-\rho)/\rho}$), 200 MC trials per cell. The seven $\rho$ levels span
weak to near-collinear within-block correlation (the grid subsumes the former second
"Rho075" preset of this demo).

---

## Output Files

Written to `simulation_results/data/`:

- `gvs_mixed_blocks_2d.txt` / `.csv` — mean FDP/TPP over the SNR $\times$ $\rho$ grid.

The file names and CSV schema match the C++ demo.

---

## Running the Demo

```bash
python demo_trex_gvs_03_mc_sim_mixed_blocks.py       # full run (hours)
python demo_trex_gvs_03_mc_sim_mixed_blocks.py 8     # with 8 worker processes
```

Monte Carlo trials are parallelized with `ProcessPoolExecutor`; the optional
first CLI argument sets the number of worker processes (default
`min(6, cpu_count)`). For a quick smoke test, temporarily reduce `num_MC` and
the sweep grids at the top of the script.

---

## Expected Results

Reference numbers and heatmaps live in the C++ counterpart's
[simulation_results](../../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_03_mc_sim_mixed_blocks/README.md#simulation-results);
the Python port reproduces the same qualitative pattern:

- **TENET** and **TENET_AUG** track each other closely and keep FDR controlled (roughly $0.05$–$0.08$, everywhere
   below target) — the inactive equicorrelated trap block does **not** inflate FDR, since its group is identified
   as null and excluded.
- **TIENET_AUG** is far more conservative (FDR near $0$, lower TPR), with a TPR collapse in the extreme
   $\rho = 0.99$ column where the within-block columns become nearly collinear.
- The random gaps and random block order mainly test that the grouping/clustering logic is invariant to block
   *placement*: TPR rises with SNR just as in the contiguous designs, so placement does not materially change
   detection power.
- Sweeping $\rho$ from $0.30$ up to $0.99$ shows how detection power and FDR control respond as within-block
   correlation strengthens, all else equal.

---

## References

1. Machkour, J., Muma, M., & Palomar, D. P., "False Discovery Rate Control for Grouped Variable Selection
   in High-Dimensional Linear Models using the T-Knock Filter.", European Signal Processing Conference (EUSIPCO), 2022,
    pp. 892–896, EURASIP.
    [DOI-Link](https://doi.org/10.23919/EUSIPCO55093.2022.9909883)
2. Machkour, J., Muma, M., & Palomar, D. P., "The Informed Elastic Net for Fast Grouped Variable Selection and
   FDR Control in Genomics Research.", Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP),
    2023, pp. 466–470, IEEE.
    [DOI-Link](https://doi.org/10.1109/CAMSAP58249.2023.10403489)
3. Zou, H., & Hastie, T. (2005). "Regularization and variable selection via the elastic net." *Journal of the Royal
   Statistical Society: Series B (Statistical Methodology)*, 67(2), pp. 301–320.
   [DOI-Link](https://doi.org/10.1111/j.1467-9868.2005.00503.x)

---

**Last updated**: 2026-07-21
