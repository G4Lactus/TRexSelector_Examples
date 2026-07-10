# Demo 01: T-Rex SPCA Monte Carlo Simulation (R)

## Purpose

Compare **T-Rex SPCA** (4 solver/mode combinations) against **ordinary PCA** and **oracle-thresholded PCA** on a synthetic sparse 3-factor model, evaluating PC1 loading-support recovery (FDR/TPR) and cumulative explained variance (PEV). This is the R counterpart of the C++ demo `demo_trex_spca_01_mc_sim.cpp` and runs the identical simulation configuration.

---

## Data Generation Parameters (`generate_sparse_factor_model()`)

$$
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E}.
$$

- $\boldsymbol{Z}$ ($n \times M$): factor columns $\mathcal{N}(0, \sigma_m^2)$ with $\sigma_m \in \{5, 3, 1\}$ for $m=1,2,3$.
- $\boldsymbol{V}$ ($p \times M$): each column has exactly `p1` nonzero entries (value $0.9$), drawn without replacement from a shared candidate pool of size `overlap_pool_size` — factor supports can partially overlap.
- $\boldsymbol{E}$: i.i.d. Gaussian, scaled to hit a target SNR in dB.
- $n=50$, $p=100$, $p_1=5$, $M=3$, `overlap_pool_size=30`.

All methods share one **center-only** $\boldsymbol{X}$ (see `center_columns()` in the suite-level `../trex_spca_sim_utils.R`: mean subtraction, no column scaling). This puts OrdPCA, OraclePCA, and TRexSPCA on a common covariance-PCA footing, matching the legacy CRAN R reference. The column scales must **not** be normalized — the factor amplitude signal lives in the column variances, and z-scoring destroys it (see the History note below).

---

## Control Parameters

```
tFDR = 0.10
lambda_2 = -1         # < 0 triggers k-fold CV for ridge-penalty selection (lambda2_method = "CV_1SE_CCD")
num_MC = 200          # single SNR point -10 dB ("TEMP fast-validation", mirroring the cpp main())
base_seed = 42        # seeds the per-trial DATA DGP deterministically:
                      #   trial mc uses base_seed_offset + mc*1000, where
                      #   base_seed_offset = base_seed + round(snr)
```

The selector itself runs with `seed = -1L` (hardware entropy for the dummies) — required for a valid Monte Carlo FDR estimate.

---

## Methods Compared

| Method | Sparsity mechanism | `mode` | `en_solver` |
|---|---|---|---|
| **OrdPCA** | None (all $p$ loadings retained) | — | — |
| **OraclePCA** | Top-$p_1$ by magnitude (true support size known) | — | — |
| **TRexSPCA-EN-Act** | T-Rex+GVS selection | `"ActiveSet"` | `"TENET"` |
| **TRexSPCA-ENaug-Act** | T-Rex+GVS selection | `"ActiveSet"` | `"TENET_AUG"` |
| **TRexSPCA-EN-Thr** | T-Rex+GVS selection | `"Thresholded"` | `"TENET"` |
| **TRexSPCA-ENaug-Thr** | T-Rex+GVS selection | `"Thresholded"` | `"TENET_AUG"` |

The PCA baselines use base-R `svd()`. Each T-Rex variant is invoked as

```r
sel <- TRexSPCASelector$new(
  X + 0, M = M, tFDR = tFDR,
  control = trex_spca_control(
    mode = "ActiveSet",          # or "Thresholded"
    gvs_type = "EN",
    en_solver = "TENET",         # or "TENET_AUG"
    lambda_2 = -1,
    lambda2_method = "CV_1SE_CCD"
  ),
  seed = -1L
)
```

with result fields `$V` (sparse loadings), `$Z` (scores), and `$active_sets` (**1-based**).

> **Binding note**: `trex_spca_control()` does not expose the per-PC selector's internal `scaling_mode`, but this is moot — the C++ demo and the Python port also use the library default (unit-L2 column scaling), which matches the legacy `lm_dummy.R` centering + L2 normalization.

---

## Metrics

- **FDR / TPR** — evaluated strictly on **PC1's** loading support only; PCs 2–3 are excluded because ordinary PCA's orthogonality constraint mixes their supports across the true factors, making a per-component ground-truth comparison ambiguous beyond PC1.
- **PEV** — cumulative percentage of explained variance across all $M$ components, via QR decomposition of the estimated score matrix.

---

## Output Files

Written to `simulation_results/` (already populated):

- `demo_trex_spca_01_mc_sim.txt` — human-readable table (methods × FDR/TPR/PEV × SNR column(s)).
- `demo_trex_spca_01_mc_sim.csv` — tidy long format: `method, metric, snr_db, value`.
- `res_demo_trex_spca_01_x.txt` — **legacy record**: recorded console output of the pre-restructure R demo (`demo_trex_spca_01.R`, now deleted), a 9-point SNR sweep with a different method set (Ord PCA, Oracle PCA, Zou `elasticnet::spca`, "Oracle SPCA", T-Rex EN/IEN Act/Thr via the CRAN-era pipeline). Since the current pipeline is back on the same centered footing, its T-Rex EN rows are roughly comparable again (e.g. FDR 10.0–17.5% at $-10$ dB) — but the method set and pipeline details differ, so treat it as historical context only.

---

## Running the Demo

```bash
Rscript R/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim.R       # default 6 worker cores
Rscript R/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim.R 8     # 8 worker cores
```

Dependencies: `TRexSelectorNeo` and base-R `parallel` only. The Monte Carlo trials run in parallel via `parallel::mclapply` (see `../trex_spca_sim_utils.R`).

---

## Real Results (SNR = -10 dB, num_MC = 200, as committed, 2026-07-08)

| Method | FDR | TPR | PEV |
|---|---|---|---|
| OrdPCA | 0.9500 | 1.0000 | 0.2017 |
| OraclePCA | 0.0040 | 0.9960 | 0.1234 |
| TRexSPCA-EN-Act | 0.1410 | 0.9990 | 0.1206 |
| TRexSPCA-ENaug-Act | 0.1288 | 0.9990 | 0.1213 |
| TRexSPCA-EN-Thr | 0.1355 | 1.0000 | 0.1328 |
| TRexSPCA-ENaug-Thr | 0.1345 | 0.9990 | 0.1331 |

Reference points: the legacy CRAN R implementation reports FDR $0.100$ / TPR $\approx 1.0$ at this same operating point. The C++ counterpart's committed table (same configuration, C++ RNG streams): OrdPCA 0.9500/1.0000/0.2021; OraclePCA 0.0060/0.9940/0.1241; T-Rex FDR $\approx 0.123$–$0.131$, TPR $\approx 0.996$–$0.997$; Python: OraclePCA 0.0080/0.9920, T-Rex FDR $\approx 0.132$–$0.136$, TPR $\approx 0.997$–$0.998$. All within Monte Carlo noise of the R table — exact numbers cannot match because the data RNG streams differ across languages. Trial-for-trial reproduction of the R table itself is also intentionally impossible: the per-trial **data** is deterministic in `base_seed`, but the dummies use hardware entropy (selector `seed = -1L`); re-runs reproduce the table to within MC noise ($\approx \pm 0.01$ at 200 trials).

> **History note (z-scoring episode, fixed 2026-07-08)**: an earlier version of this pipeline z-score-standardized $\boldsymbol{X}$, which converts the covariance PCA into a correlation PCA and destroys the factor amplitude signal — the measured T-Rex FDR inflated to $\approx 0.53$ and OraclePCA TPR dropped to $\approx 0.84$ at this same nominal SNR. The C++ validation program `validation_trex_spca_06_handrolled_comparison` isolated the cause (the `TRexSPCA` class and a hand-rolled pipeline replicating the legacy R recipe agree within MC noise; only the preprocessing differed), and the pipeline is back to center-only.

---

## Interpretation

- **OrdPCA**'s $\mathrm{FDR}\approx0.95$ and $\mathrm{TPR}=1.0$ are trivial artifacts of retaining all $p=100$ loadings (true support is $p_1=5$) — not a genuine comparison point, only a non-sparse reference.
- **OraclePCA** achieves near-perfect PC1 recovery (it is told the true support size); its PEV ($0.123$) is the useful reference value here — how much variance the true sparse signal explains, independent of any estimation procedure.
- All four **T-Rex SPCA** variants recover PC1's support almost perfectly (TPR $\approx 1.0$) with FDR close to the $\mathrm{tFDR}=0.10$ target — even though $-10$ dB sounds like a very weak signal, the strong first factor (std $5$) dominates PC1 recovery on the covariance-PCA footing.
- The realized FDR ($\approx 0.13$–$0.14$) sits $\approx 0.03$ above the target and the legacy CRAN R reference ($0.100$). This is the long-known residual gap between the C++ core and the CRAN R implementation, investigated by the programs now in the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`) (suspects: the internal CD-ridge λ₂ CV vs. R's `cv.glmnet`, and dummy-generation details); `validation_trex_spca_06` confirms it is **not** a `TRexSPCASelector` artifact.
- The four T-Rex variants are statistically indistinguishable from each other at this single SNR point; no solver (`TENET` vs. `TENET_AUG`) or mode (`ActiveSet` vs. `Thresholded`) conclusion should be drawn from it. The **Thresholded** variants' slightly higher PEV ($\approx 0.133$ vs. $\approx 0.121$) is the one systematic difference: they rebuild loadings from the full ordinary-PCA solution on the selected support, which preserves more explained variance than the active-set coefficients.

---

**Last updated**: 2026-07-08
