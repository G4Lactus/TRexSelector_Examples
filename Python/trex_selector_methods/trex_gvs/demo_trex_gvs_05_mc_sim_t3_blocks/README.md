# Demo 05: T-Rex+GVS on Heavy-Tailed (Student-t(3)) Blocks (Python)

## Purpose

Stress-test T-Rex+GVS on non-Gaussian data.
 This demo reuses the same four-block active/trap geometry as Demo 03 (Mixed Blocks), but replaces every Gaussian
 component of the data-generating process — latent factors, within-block perturbations, and response noise — with
 scaled Student-t(3) variables, testing robustness to heavy-tailed departures from the Gaussian assumption.
 Selection and evaluation are per variable (see
 [What is actually measured](../README.md#what-is-actually-measured-in-these-demos)).

This is the Python port of the C++ demo
[demo_trex_gvs_05_mc_sim_t3_blocks](../../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_05_mc_sim_t3_blocks/),
using `trex_selector_neo`'s `TRexGVSSelector`.

---

## Data Generation Parameters (`dgp_t3_equi_blocks_snr`)

We consider again the linear model, but with heavy-tailed rather than Gaussian noise:

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \sigma_{\varepsilon}\, \boldsymbol{\epsilon},
\qquad \epsilon_i \sim t(3)/\sqrt{3}
$$

- $\boldsymbol{y} \in \mathbb{R}^n$ is the response vector.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ is the design matrix.
- $\boldsymbol{\beta} \in \mathbb{R}^p$ is the coefficient vector, with $s$ nonzero entries.
- $\boldsymbol{\epsilon}$ is the noise vector, i.i.d. $t(3)/\sqrt{3}$ — a Student-t distribution with 3 degrees of
   freedom, rescaled to unit variance.
- $\sigma_{\varepsilon}^2$ is the noise variance, calibrated to achieve a target linear signal-to-noise ratio (SNR).
- $n = 200$, $p = 500$, $s = 150$ (the sum of the three active blocks).

Unlike Demos 01–04, the design matrix uses the **unit-variance** equicorrelated form, parameterized directly by
 $\rho$ rather than through $\sigma_x$:

$$
X_{ij} = \sqrt{\rho}\, Z_{i,\,g(j)} + \sqrt{1-\rho}\, \xi_{ij},
\qquad Z_{\cdot,k},\, \xi_{ij} \sim t(3)/\sqrt{3},
$$

- Same 4-block geometry as Demo 03: sizes $\{20, 50, 80, 65\}$ (215 grouped columns in total).
- **Blocks 0–2 are active** ($\beta_j = 3$); **block 3 is an inactive heavy-tailed trap** ($\beta_j = 0$).
- Both the latent factors $Z_{\cdot,k}$ and the within-block perturbations $\xi_{ij}$ are $t(3)/\sqrt{3}$, so the
   heavy tails enter the design matrix as well as the response.

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
$\rho \in \{0.30, 0.50, 0.70, 0.85, 0.90, 0.99\}$, 200 MC trials per cell. Here $\rho$ enters the
DGP directly (no $\sigma_x$ reparameterization).

---

## Output Files

Written to `simulation_results/data/`:

- `gvs_t3_blocks_2d.txt` / `.csv` — mean FDP/TPP over the SNR $\times$ $\rho$ grid.

The file names and CSV schema match the C++ demo.

---

## Running the Demo

```bash
python demo_trex_gvs_05_mc_sim_t3_blocks.py       # full run (hours)
python demo_trex_gvs_05_mc_sim_t3_blocks.py 8     # with 8 worker processes
```

Monte Carlo trials are parallelized with `ProcessPoolExecutor`; the optional
first CLI argument sets the number of worker processes (default
`min(6, cpu_count)`). For a quick smoke test, temporarily reduce `num_MC` and
the sweep grids at the top of the script.

---

## Expected Results

Reference numbers and heatmaps live in the C++ counterpart's
[simulation_results](../../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_05_mc_sim_t3_blocks/README.md#simulation-results);
the Python port reproduces the same qualitative pattern:

- Despite the $t(3)$ noise, FDR control holds throughout: **TENET** and **TENET_AUG** (identical here) stay just
   under the $\mathrm{tFDR} = 0.1$ target in every cell, with TPR rising steadily with SNR.
- The heavy tails cost the TENET variants a little power relative to the equivalent Gaussian case (Demo 03) at
   comparable SNR — occasional large outliers make exact group boundaries noisier to detect — but they do not
   compromise FDR control. That is the headline result: heavy tails buy lower TPR, not FDR drift.
- **TIENET_AUG** is the most conservative of the three within this demo, and its power collapses entirely in the
   extreme $\rho = 0.99$ column, where the within-block columns become nearly collinear.
- Judged against *its own* behavior on the Gaussian equicorrelated designs (Demos 01–03), though, TIENET_AUG is far
   more competitive here: away from the $\rho = 0.99$ column its TPR is substantial rather than near-zero.

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
