# T-Rex SPCA: Sparse Principal Component Analysis Demonstration Suite

## Overview

This folder contains C++ example programs for **T-Rex SPCA**, a sparse principal component analysis method
 built on the T-Rex+GVS elastic-net machinery (see [../trex_gvs/](../trex_gvs/README.md)) according to
 [[1,2]](#references).

Instead of selecting sparse *variables* in a regression model, T-Rex SPCA selects sparse *loadings* for each
 principal component: the T-Rex+GVS selector is applied per component to control the false discovery rate of
 the estimated loading support. The result is a PCA whose components are supported on a small,
 FDR-controlled set of variables rather than on all $p$ of them.

The demo in this folder is designed to help users understand:

1. whether the FDR of the estimated loading support is actually held at the target,
2. what that control costs in explained variance, measured against ordinary PCA and an oracle that is told
    the true support size,
3. how the elastic-net solver (`TENET` vs. `TENET_AUG`) and the loading-assembly mode (`ActiveSet` vs.
    `Thresholded`) affect the outcome.

---

## What this folder covers

- A **sparse $M$-factor data-generating process** where each column of the loading matrix has exactly $p_1$
   nonzero entries drawn from a shared overlap pool, so factor supports can partially coincide.
- Two **T-Rex SPCA** solver variants crossed with two loading-assembly modes, benchmarked against **ordinary
   PCA** (no sparsity) and **oracle-thresholded PCA** (knows the true support size).
- A sweep over the signal-to-noise ratio **in decibels**, the scale this suite reports on throughout.

---

## Statistical model and targets

### The factor model

```math
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E},
\qquad
\boldsymbol{Z}_{\cdot,m} \sim \mathcal{N}(\boldsymbol{0}, \sigma_m^2 \boldsymbol{I}_n).
```

### Notation

- $n$: number of observations, $p$: number of variables, $M$: number of latent factors.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$: centered observation matrix.
- $\boldsymbol{V} \in \mathbb{R}^{p \times M}$: sparse loading matrix; each column carries exactly $p_1$
   nonzero entries of value $0.9$, drawn without replacement from a shared candidate pool of size
   `overlap_pool_size`.
- $\sigma_m$: amplitude of factor $m$, decreasing over $\{5, 3, 1\}$ (falling back to $1$ for $M > 3$), so
   the components are ordered by strength.
- $\mathcal{S}_1$: PC1's true active loading support; $\widehat{\mathcal{S}}_1$ its estimated counterpart.
- $\boldsymbol{E}$: i.i.d. Gaussian noise, scaled to hit the target SNR.

### Signal-to-noise ratio in decibels

```math
\mathrm{SNR}_{\mathrm{dB}} = 10 \log_{10}\!\left(\frac{\mathrm{Var}(\boldsymbol{Z}\boldsymbol{V}^\top)}{\mathrm{Var}(\boldsymbol{E})}\right).
```

Both the tables and the figures use this decibel axis directly; it is never converted back to a linear ratio.

### False discovery rate, true positive rate, explained variance

```math
\mathrm{FDR} = \mathbb{E}\left[\frac{|\widehat{\mathcal{S}}_1 \setminus \mathcal{S}_1|}{\max\{1, |\widehat{\mathcal{S}}_1|\}}\right],
\qquad
\mathrm{TPR} = \mathbb{E}\left[\frac{|\mathcal{S}_1 \cap \widehat{\mathcal{S}}_1|}{\max\{1, |\mathcal{S}_1|\}}\right].
```

**PEV** — the proportion of explained variance, computed cumulatively across all $M$ components via a QR
 decomposition of the estimated score matrix — is the quantity the sparsity is bought against. A method that
 selects fewer loadings explains less variance, so PEV is where the cost of FDR control becomes visible.

### What is actually measured in these demos

> **FDR and TPR are evaluated on the first principal component only.** Components 2 and 3 are excluded
> deliberately: ordinary PCA's orthogonality constraint mixes their loading supports across the true factors,
> so there is no unambiguous per-component ground truth beyond PC1. This is a scope limitation of the
> evaluation, not of T-Rex SPCA. PEV does not suffer the same ambiguity and is therefore cumulative over all
> $M$ components.
>
> Note also that **ordinary PCA trivially attains $\mathrm{TPR} = 1$ and $\mathrm{FDR} \approx (p - p_1)/p$**
> by retaining all $p$ loadings — a non-sparse reference point, not a competitor on these metrics. **Oracle
> PCA**, told the true support cardinality $p_1$, bounds what any data-driven sparse method could achieve.
>
> Every method is evaluated on the same **center-only** design matrix — mean subtraction, no column scaling —
> which puts them on a common covariance-PCA footing. Column scaling must be avoided here: the factor
> amplitudes live in the column variances, so z-scoring (correlation PCA) destroys the very signal that
> distinguishes the factors.

---

## T-Rex SPCA specific concepts

- **`SPCAMode`** — how the sparse loading vector is assembled once a support has been selected:
  - `ActiveSet`: use only the T-Rex+GVS selected variables' loadings,
  - `Thresholded`: assemble loadings from the full ordinary-PCA solution, restricted to the selected support.
- **`ENSolverType`** — the elastic-net solver driving each per-component selection: `TENET` (Gram-based) or
   `TENET_AUG` (augmented-LASSO formulation); see [../trex_gvs/](../trex_gvs/README.md) for the shared
   elastic-net background.
- **`lambda_2`** — the ridge penalty: a negative value selects it by $k$-fold cross-validation, $0$ disables
   the ridge term, and a positive value fixes it.
- **`tFDR`** — the target FDR handed to the per-component T-Rex selector.

---

## Folder contents

```txt
trex_spca/
  ├── README.md
  ├── CMakeLists.txt
  ├── trex_spca_sim_utils.hpp          # DGP, MC loop, table + CSV output
  ├── trex_spca_plt_utils.py           # shared plotting module (dB sweep figures)
  │
  └── demo_trex_spca_01_mc_sim/
      ├── demo_trex_spca_01_mc_sim.cpp # the demo source
      ├── README.md                    # scenario description and results
      ├── generate_plots.sh            # renders this demo's figure from its CSV
      └── simulation_results/
          ├── data/                    # .txt summary + .csv table
          └── plots/                   # .png + .pdf figures
```

Correctness probes for the individual pipeline stages — TENET vs. TENET_AUG equivalence, scaling
 comparisons, step-by-step C++/R diffs, and a dump of the GVS selector's internal calibration matrices —
 live in the library's test suite under
 `TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`, built via its `build_diagnostics`
 target. The R and Python counterparts of this demo run the identical simulation and can be used to
 cross-check the numbers.

---

## Building the demo

From the C++ workspace root:

```bash
cd TRexSelector_Examples/cpp
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"
cmake --build build/release
```

The executable mirrors the source tree under the `bin/` directory:

```txt
build/release/bin/trex_selector_methods/trex_spca/<demo_folder>/<demo_name>
```

---

## Running a demo

```bash
./build/release/bin/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim
```

The demo writes a human-readable `.txt` summary and a tidy `.csv` table into its local
`simulation_results/data/` folder.

Note that only the data is seeded deterministically; each trial's dummies are drawn from hardware entropy by
 design, which is what makes the Monte Carlo FDR estimate valid. A re-run therefore reproduces the reported
 behaviour to within Monte Carlo noise rather than exactly.

### Figures

The demo folder has a `generate_plots.sh` that renders its figure from the saved CSV into
`simulation_results/plots/` using the shared `trex_spca_plt_utils.py` and the repo-root `.venv`
(matplotlib + pandas + numpy):

```bash
cd demo_trex_spca_01_mc_sim
./generate_plots.sh
```

The figure shows TPR, FDR (against the tFDR target) and PEV over the decibel sweep axis. `--xscale` is
accepted for parity with the other suites, but the axis stays in dB either way.

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

---

**Last updated**: 2026-07-20
