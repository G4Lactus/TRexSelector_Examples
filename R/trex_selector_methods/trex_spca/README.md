# T-Rex SPCA: Sparse Principal Component Analysis — R Demonstration Suite

## Overview

This folder contains the R example program for **T-Rex SPCA**, a sparse principal component analysis method built on top of **T-Rex+GVS** (see [../trex_gvs/](../trex_gvs/README.md)). Instead of selecting sparse *variables* in a regression model, T-Rex SPCA selects sparse *loadings* for each principal component, applying the T-Rex+GVS elastic-net machinery per component to control the false discovery rate of the estimated loading support.

The single demo in this folder compares T-Rex SPCA against two PCA baselines on a synthetic sparse-factor-model dataset. It is the R counterpart of the C++ demo `demo_trex_spca_01_mc_sim.cpp` and mirrors that suite's folder structure one-to-one.

---

## What this folder covers

- A **sparse M-factor data-generating process** $X = ZV^\top + E$, where the loading matrix $V$ has exactly `p1` nonzero entries per factor, drawn from a shared overlap pool so factor supports can partially coincide.
- Two **T-Rex SPCA** solver/mode combinations (elastic-net solver `TENET` vs. `TENET_AUG`; loading-assembly mode `ActiveSet` vs. `Thresholded`) benchmarked against **ordinary PCA** (no sparsity) and **oracle-thresholded PCA** (knows the true support size). The PCA baselines use base-R `svd()`.

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

All methods are evaluated on the same $\boldsymbol{X}$ after a shared **center-only** preprocessing step (mean subtraction, no column scaling — a covariance-PCA footing), so that OrdPCA, OraclePCA, and T-Rex SPCA are compared on equal preprocessing. The column scales must **not** be normalized: in this factor model the amplitude signal lives in the column variances, and z-scoring the columns (correlation PCA) destroys it — measured effect at $-10$ dB: T-Rex SPCA FDR degrades from $\approx 0.13$ to $\approx 0.53$ (see the C++ `validation_trex_spca_06_handrolled_comparison`).

---

## T-Rex SPCA specific concepts

- **`mode`** (in `trex_spca_control()`) — how the sparse loading vector for each principal component is assembled from the underlying T-Rex+GVS selection:
  - `"ActiveSet"`: use only the T-Rex+GVS *selected* variables' loadings.
  - `"Thresholded"`: assemble loadings from the full ordinary-PCA solution, thresholded down to the T-Rex+GVS-selected support.
- **`en_solver`** — the elastic-net solver driving each per-component T-Rex+GVS selection: `"TENET"` (Gram-based) or `"TENET_AUG"` (augmented-LASSO formulation); see [../trex_gvs/README.md](../trex_gvs/README.md) for the shared elastic-net background.
- **PC1-only FDR/TPR** — FDR and TPR are evaluated **strictly on the first principal component's loading support**. PCs 2 and 3 are intentionally excluded from these metrics because ordinary PCA's orthogonality constraint mixes their loading supports across the true factors, making a per-component "true support" comparison ambiguous beyond PC1. This is a deliberate scope limitation of the evaluation, not a limitation of T-Rex SPCA itself.
- **Cumulative PEV** — percentage of explained variance is instead evaluated cumulatively across **all** $M$ components (via a QR decomposition of the estimated score matrix), since this metric does not suffer from the same per-component ambiguity.
- **Internal scaling (R binding note)** — `trex_spca_control()` does not expose the per-PC selector's internal `scaling_mode`, but this is moot: the C++ demo and the Python port also use the library default (unit-L2 column scaling), which matches the legacy `lm_dummy.R` centering + L2 normalization.

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

where $\mathcal{S}_1$ is PC1's true active loading support and $\widehat{\mathcal{S}}_1$ its estimated counterpart (both **1-based** in R). **Ordinary PCA** trivially achieves $\mathrm{TPR}=1$ and $\mathrm{FDR}\approx(p-p_1)/p$ since it retains all $p$ loadings — it is a non-sparse reference point, not a genuine competitor on these two metrics. **Oracle PCA** (knows the true support cardinality $p_1$) provides an upper bound on what any data-driven sparse method could achieve.

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
  ├── trex_spca_sim_utils.R
  └── demo_trex_spca_01_mc_sim/
      ├── demo_trex_spca_01_mc_sim.R
      ├── README.md
      └── simulation_results/
          ├── demo_trex_spca_01_mc_sim.txt
          ├── demo_trex_spca_01_mc_sim.csv
          └── res_demo_trex_spca_01_x.txt   <- legacy pre-restructure record
```

The suite-level [trex_spca_sim_utils.R](trex_spca_sim_utils.R) mirrors the C++ `trex_spca_sim_utils.hpp`: `generate_sparse_factor_model()`, `center_columns()`, `evaluate_spca()`, the Monte Carlo runners `run_mc_trials_pca()` / `run_mc_trials_oracle_pca()` / `run_mc_trials_trex_spca()` (parallel via `parallel::mclapply`), and `save_and_print_spca_mc_results()` (console + `.txt` table plus a tidy CSV `method,metric,snr_db,value`).

---

## Example results

The demo has **real committed output** — a single SNR point ($-10$ dB), `num_MC = 200`, produced by the center-only pipeline:

| Method | FDR | TPR | PEV |
|---|---|---|---|
| OrdPCA | 0.9500 | 1.0000 | 0.2017 |
| OraclePCA | 0.0040 | 0.9960 | 0.1234 |
| TRexSPCA-EN-Act | 0.1410 | 0.9990 | 0.1206 |
| TRexSPCA-ENaug-Act | 0.1288 | 0.9990 | 0.1213 |
| TRexSPCA-EN-Thr | 0.1355 | 1.0000 | 0.1328 |
| TRexSPCA-ENaug-Thr | 0.1345 | 0.9990 | 0.1331 |

### Interpretation

- **OrdPCA**'s FDR $\approx 0.95$ is the expected trivial artifact of selecting all $p=100$ loadings (true support is only $p_1=5$) — not a meaningful comparison point.
- **OraclePCA** achieves near-perfect recovery (it is told the true support size); its PEV ($0.123$) is the reference for how much variance the true sparse signal explains, independent of any selection procedure.
- All four **T-Rex SPCA** variants recover PC1's support almost perfectly (TPR $\approx 1.0$) with FDR close to the $\mathrm{tFDR}=0.10$ target. The realized FDR ($\approx 0.13$–$0.14$) sits $\approx 0.03$ above the target and the legacy CRAN R reference ($0.100$ at this point) — the long-known residual gap between the C++ core and the CRAN R implementation, investigated by the cross-check programs in the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`); the C++ `validation_trex_spca_06_handrolled_comparison` confirms it is not a `TRexSPCASelector` artifact.
- The four T-Rex variants are statistically indistinguishable from each other at this one data point; no solver/mode conclusion should be drawn from a single SNR level. The Thresholded variants' slightly higher PEV is the one systematic difference.
- The C++ counterpart's committed table (same configuration, C++ RNG streams) agrees within Monte Carlo noise: OrdPCA 0.9500/1.0000/0.2021, OraclePCA 0.0060/0.9940/0.1241, T-Rex FDR $\approx 0.123$–$0.131$, TPR $\approx 0.996$–$0.997$; likewise Python (T-Rex FDR $\approx 0.132$–$0.136$). Exact numbers cannot match across languages because the data RNG streams differ.

> **History note (z-scoring episode, fixed 2026-07-08)**: an earlier version of this pipeline z-score-standardized $\boldsymbol{X}$, which destroyed the factor amplitude signal and inflated the measured T-Rex FDR to $\approx 0.53$ at this same nominal SNR; the C++ `validation_trex_spca_06_handrolled_comparison` isolated the cause and the pipeline is back to center-only. Like the C++ `main()`, the R demo currently exercises a single SNR value ($-10$ dB, "TEMP fast-validation") with `num_MC = 200`. Trial-for-trial reproduction is intentionally impossible: each trial's **data** is drawn from a deterministic seed, but the per-trial dummies use hardware entropy (selector `seed = -1L`), as required for a valid Monte Carlo FDR estimate; re-runs reproduce the table to within MC noise. See the demo's own [README.md](demo_trex_spca_01_mc_sim/README.md) for detail.

---

## Running the demo

```bash
Rscript R/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim.R       # default 6 worker cores
Rscript R/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim.R 8     # 8 worker cores
```

Dependencies: `TRexSelectorNeo` and base-R `parallel` only (the old `elasticnet`/`glmnet` dependencies are gone from this suite).

---

## Cross-check probes (relocated)

The R↔C++ cross-check material that used to live in this folder — `demo_trex_spca_02.R` (R reference-dump generator), `lambda2_probe.R`, `lambda2_foldmatch.R`, `probe_R_dummy_variance.R`, `r_dummy_variance_probe.R`, and the baked `rdump/` / `rdump10/` data folders — moved to the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`). The `rdump_snr7_backup/` folder was intentionally dropped, not carried over. See that suite for the stage-by-stage pipeline comparisons; this README no longer reproduces the old probe tables.

---

## Notes for new users

- Start with `demo_trex_spca_01_mc_sim/` — it is currently the only demo in this folder.
- Read the PC1-only FDR/TPR rationale above before interpreting results — this is a deliberate evaluation-scope choice, not a bug.
- Support sets and `$active_sets` indices are **1-based** in R (the C++ sources are 0-based).
- Counterparts: the C++ suite at [../../../cpp/trex_selector_methods/trex_spca/](../../../cpp/trex_selector_methods/trex_spca/); the Python suite at [../../../Python/trex_selector_methods/trex_spca/](../../../Python/trex_selector_methods/trex_spca/); and the C++/R validation probes, now in the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`).

---

## References

- See [../README.md](../README.md) for the category overview of all T-Rex selector variants.
- See [../trex_gvs/README.md](../trex_gvs/README.md) for the shared T-Rex+GVS elastic-net background (`gvs_type`, `en_solver`, `lambda2_method`).

---

**Last updated**: 2026-07-08
