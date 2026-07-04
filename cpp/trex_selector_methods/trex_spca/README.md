# T-Rex SPCA: Sparse Principal Component Analysis Demonstration Suite

## Overview

This folder contains C++ example programs for **T-Rex SPCA**, a sparse principal component analysis method built on top of **T-Rex+GVS** (see [../trex_gvs/](../trex_gvs/README.md)). Instead of selecting sparse *variables* in a regression model, T-Rex SPCA selects sparse *loadings* for each principal component, applying the T-Rex+GVS elastic-net machinery per component to control the false discovery rate of the estimated loading support.

The single demo in this folder compares T-Rex SPCA against two PCA baselines on a synthetic sparse-factor-model dataset.

---

## What this folder covers

- A **sparse M-factor data-generating process** $X = ZV^\top + E$, where the loading matrix $V$ has exactly `p1` nonzero entries per factor, drawn from a shared overlap pool so factor supports can partially coincide.
- Two **T-Rex SPCA** solver/mode combinations (elastic-net solver `TENET` vs. `TENET_AUG`; loading-assembly mode `ActiveSet` vs. `Thresholded`) benchmarked against **ordinary PCA** (no sparsity) and **oracle-thresholded PCA** (knows the true support size).

---

## Statistical model and notation

$$
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E},
$$

- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$: centered observation matrix.
- $\boldsymbol{Z} \in \mathbb{R}^{n \times M}$: true latent factor scores, column $m$ drawn $\mathcal{N}(0, \sigma_m^2)$ with decreasing factor strength $\sigma_m \in \{5, 3, 1\}$ (falling back to $1$ for $M > 3$).
- $\boldsymbol{V} \in \mathbb{R}^{p \times M}$: true sparse loading matrix; each column has exactly `p1` nonzero entries (value $0.9$), drawn without replacement from a shared candidate pool of size `overlap_pool_size` — so different factors' active sets can overlap.
- $\boldsymbol{E} \in \mathbb{R}^{n \times p}$: i.i.d. Gaussian noise, scaled so that the realized signal-to-noise ratio matches a target value **in decibels**,
  $$
  \mathrm{SNR}_{\mathrm{dB}} = 10 \log_{10}\!\left(\frac{\mathrm{Var}(\boldsymbol{Z}\boldsymbol{V}^\top)}{\mathrm{Var}(\boldsymbol{E})}\right).
  $$

All methods are evaluated on the same $\boldsymbol{X}$ after a shared z-score column standardization (a correlation-PCA footing), so that OrdPCA, OraclePCA, and T-Rex SPCA are compared on equal preprocessing.

---

## T-Rex SPCA specific concepts

- **`SPCAMode`** — how the sparse loading vector for each principal component is assembled from the underlying T-Rex+GVS selection:
  - `ActiveSet`: use only the T-Rex+GVS *selected* variables' loadings.
  - `Thresholded`: assemble loadings from the full ordinary-PCA solution, thresholded down to the T-Rex+GVS-selected support.
- **`ENSolverType`** — the elastic-net solver driving each per-component T-Rex+GVS selection: `TENET` (Gram-based) or `TENET_AUG` (augmented-LASSO formulation); see [../trex_gvs/README.md](../trex_gvs/README.md) for the shared elastic-net background.
- **PC1-only FDR/TPR** — FDR and TPR are evaluated **strictly on the first principal component's loading support**. PCs 2 and 3 are intentionally excluded from these metrics because ordinary PCA's orthogonality constraint mixes their loading supports across the true factors, making a per-component "true support" comparison ambiguous beyond PC1. This is a deliberate scope limitation of the evaluation, not a limitation of T-Rex SPCA itself.
- **Cumulative PEV** — percentage of explained variance is instead evaluated cumulatively across **all** $M$ components (via a QR decomposition of the estimated score matrix), since this metric does not suffer from the same per-component ambiguity.

---

## Statistical targets

$$
\mathrm{FDR}
=
\mathbb{E}
\left[
\frac{|\widehat{\mathcal{S}}_1 \setminus \mathcal{S}_1|}
{\max\{1, |\widehat{\mathcal{S}}_1|\}}
\right],
\qquad
\mathrm{TPR}
=
\mathbb{E}
\left[
\frac{|\mathcal{S}_1 \cap \widehat{\mathcal{S}}_1|}
{\max\{1, |\mathcal{S}_1|\}}
\right],
$$

where $\mathcal{S}_1$ is PC1's true active loading support and $\widehat{\mathcal{S}}_1$ its estimated counterpart. **Ordinary PCA** trivially achieves $\mathrm{TPR}=1$ and $\mathrm{FDR}\approx(p-p_1)/p$ since it retains all $p$ loadings — it is a non-sparse reference point, not a genuine competitor on these two metrics. **Oracle PCA** (knows the true support cardinality $p_1$) provides an upper bound on what any data-driven sparse method could achieve.

---

## Start here

There is currently one demo: **`demo_trex_spca_01_mc_sim/`**, comparing all methods via a Monte Carlo simulation at a fixed SNR. See its local [README.md](demo_trex_spca_01_mc_sim/README.md) for the full parameter breakdown and the real committed result numbers.

---

## Demo descriptions

| # | Name | Purpose | Setting | Main parameters |
|---|------|---------|---------|-----------------|
| **01** | MC Sim | Compare T-Rex SPCA (2 solvers × 2 modes) against OrdPCA/OraclePCA | Sparse 3-factor model, $n=50$, $p=100$ | $p_1=5$, `overlap_pool_size=30`, tFDR=0.10, currently run at a single SNR point |

---

## Folder contents

```txt
trex_spca/
  ├── README.md
  ├── CMakeLists.txt
  ├── trex_spca_sim_utils.hpp
  └── demo_trex_spca_01_mc_sim/
      ├── demo_trex_spca_01_mc_sim.cpp
      ├── README.md
      └── simulation_results/
          ├── demo_trex_spca_01_mc_sim.txt
          └── demo_trex_spca_01_mc_sim.csv
```

---

## Example results

Unlike some of the other T-Rex selector variant folders, this demo already has **real committed output** — from a single SNR point ($-10$ dB):

| Method | FDR | TPR | PEV |
|---|---|---|---|
| OrdPCA | 0.9500 | 1.0000 | 0.2024 |
| OraclePCA | 0.0000 | 1.0000 | 0.1254 |
| TRexSPCA-EN-Act | 0.1188 | 0.9975 | 0.1204 |
| TRexSPCA-ENaug-Act | 0.1337 | 0.9975 | 0.1205 |
| TRexSPCA-EN-Thr | 0.1099 | 1.0000 | 0.1252 |
| TRexSPCA-ENaug-Thr | 0.1107 | 1.0000 | 0.1253 |

### Interpretation

- **OrdPCA**'s FDR $\approx 0.95$ is the expected trivial artifact of selecting all $p=100$ loadings (true support is only $p_1=5$) — not a meaningful comparison point.
- **OraclePCA** achieves perfect recovery by construction (it is told the true support size); its PEV ($0.1254$) is a useful reference for "how much variance the true sparse signal explains" independent of any selection procedure.
- All four **T-Rex SPCA** variants recover PC1's support almost perfectly (TPR $\approx 0.998$–$1.0$) while keeping FDR close to the $\mathrm{tFDR}=0.10$ target — even though $-10$ dB sounds like a very weak signal, the shared-factor construction (factor scale $5$ for the first factor) makes PC1 easy to recover in this particular design.
- **Thresholded** mode is very slightly ahead of **ActiveSet** mode at this one data point (marginally higher TPR, FDR marginally closer to target) — this is a single-point observation from one SNR level and should not be read as a general conclusion about the two modes.

> **Caveat**: the demo's own file-level doc-comment describes a 5-point SNR sweep ($\{-10,-5,0,5,10\}$ dB, num_MC=100), but the `main()` function currently only exercises a single SNR value ($-10$ dB, explicitly commented `// TEMP fast-validation`) with `num_MC=200` — while the committed output file header reports "averaged over 80 MC trials". These three numbers do not agree with each other; treat the table above as illustrative of this one SNR point under whatever configuration actually produced the committed file, not as a reproducible-on-rebuild guarantee. See the demo's own [README.md](demo_trex_spca_01_mc_sim/README.md) for detail.

---

## Building the demo

```bash
cd TRexSelector_Simulations/cpp
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"
cmake --build build/debug
```

---

## Running the demo

```bash
./build/debug/bin/demo_trex_spca_01_mc_sim
```

---

## Notes for new users

- Start with `demo_trex_spca_01_mc_sim/` — it is currently the only demo in this folder.
- Read the PC1-only FDR/TPR rationale above before interpreting results — this is a deliberate evaluation-scope choice, not a bug.
- Cross-check numbers against the R reference implementation in `R/trex_selector_methods/trex_spca/` (`demo_trex_spca_01.R`, `demo_trex_spca_02.R`, plus `lambda2_probe.R` / `lambda2_foldmatch.R` for the ridge-penalty selection step specifically).
- For deeper correctness checks of individual pipeline stages, see the five validation programs in `cpp/trex_selector_methods/validation/trex_spca/`: `validation_trex_spca_01_lambda2_probe.cpp` (λ₂ CV vs. R `cv.glmnet`), `validation_trex_spca_02_solver_comparison.cpp` (TENET vs. TENET_AUG equivalence), `validation_trex_spca_03_scaling_comparison.cpp` (unit-L2 vs. z-score scaling), `validation_trex_spca_04_rdump_pipeline.cpp` (step-by-step C++/R diff), and `validation_trex_spca_05_lambda2_foldmatch.cpp` (identical-fold CV curve comparison).

---

## References

- See [../README.md](../README.md) for the category overview of all T-Rex selector variants.
- See [../trex_gvs/README.md](../trex_gvs/README.md) for the shared T-Rex+GVS elastic-net background (`ENSolverType`, `LambdaSelectionMethod`).

---

**Last updated**: 2026-07-04
