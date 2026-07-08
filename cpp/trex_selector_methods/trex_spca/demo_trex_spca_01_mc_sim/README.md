# Demo 01: T-Rex SPCA Monte Carlo Simulation

## Purpose

Compare **T-Rex SPCA** (4 solver/mode combinations) against **ordinary PCA** and **oracle-thresholded PCA** on a synthetic sparse 3-factor model, evaluating PC1 loading-support recovery (FDR/TPR) and cumulative explained variance (PEV).

---

## Data Generation Parameters (`DataGenerator::generate_sparse_factor_model`)

$$
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E}.
$$

- $\boldsymbol{Z}$ ($n \times M$): factor columns $\mathcal{N}(0, \sigma_m^2)$ with $\sigma_m \in \{5, 3, 1\}$ for $m=0,1,2$.
- $\boldsymbol{V}$ ($p \times M$): each column has exactly `p1` nonzero entries (value $0.9$), drawn without replacement from a shared candidate pool of size `overlap_pool_size` — factor supports can partially overlap.
- $\boldsymbol{E}$: i.i.d. Gaussian, scaled to hit a target SNR in dB.
- $n=50$, $p=100$, $p_1=5$, $M=3$, `overlap_pool_size=30`.

All methods share one **center-only** $\boldsymbol{X}$ (see `center_columns()`: mean subtraction, no column scaling). This puts OrdPCA, OraclePCA, and TRexSPCA on a common covariance-PCA footing, matching the legacy R reference. The column scales must **not** be normalized — the factor amplitude signal lives in the column variances, and z-scoring destroys it (see the History note below).

---

## Control Parameters

```
tFDR = 0.10
lambda_2 = -1        # < 0 triggers k-fold CV for ridge-penalty selection (0 = no ridge, > 0 = fixed)
num_MC = 200          # as set in main(), see caveat below
base_seed = 42        # seeds the per-trial DATA DGP deterministically (trial mc uses base_seed + mc*1000)
```

---

## Methods Compared

| Method | Sparsity mechanism | `SPCAMode` | `ENSolverType` |
|---|---|---|---|
| **OrdPCA** | None (all $p$ loadings retained) | — | — |
| **OraclePCA** | Top-$p_1$ by magnitude (true support size known) | — | — |
| **TRexSPCA-EN-Act** | T-Rex+GVS selection | `ActiveSet` | `TENET` |
| **TRexSPCA-ENaug-Act** | T-Rex+GVS selection | `ActiveSet` | `TENET_AUG` |
| **TRexSPCA-EN-Thr** | T-Rex+GVS selection | `Thresholded` | `TENET` |
| **TRexSPCA-ENaug-Thr** | T-Rex+GVS selection | `Thresholded` | `TENET_AUG` |

All four T-Rex SPCA variants run their per-component selector with the default `ScalingMode::L2` (matching the legacy `lm_dummy.R` column centering + L2 normalization).

---

## Metrics

- **FDR / TPR** — evaluated strictly on **PC1's** loading support only; PCs 2–3 are excluded because ordinary PCA's orthogonality constraint mixes their supports across the true factors, making a per-component ground-truth comparison ambiguous beyond PC1.
- **PEV** — cumulative percentage of explained variance across all $M$ components, via QR decomposition of the estimated score matrix.

---

## ⚠️ Configuration notes

`main()` currently sets `cfg.num_MC = 200` and `snr_values = {-10.0}` — a **single SNR point**, explicitly commented `// TEMP fast-validation`. Only this one SNR point has been run and committed; the file's doc-comment now describes exactly this configuration.

The committed output (`simulation_results/`) matches the current `main()` configuration (200 trials, $-10$ dB) and the center-only pipeline. **History note**: during an earlier debugging round the pipeline briefly gained a shared z-score standardization; that converts the covariance PCA into a correlation PCA, destroys the factor amplitude signal, and pushed the measured T-Rex FDR to $\approx 0.52$ (OraclePCA TPR to $\approx 0.87$) at this same nominal SNR. `validation_trex_spca_06_handrolled_comparison` isolated the cause (the class and a hand-rolled legacy-recipe pipeline agree; only the preprocessing differed) and the z-scoring was removed on 2026-07-08.

Trial-for-trial reproduction is intentionally impossible: each trial's **data** is drawn from a deterministic seed (`base_seed + mc*1000`), so the sequence of design matrices is fixed by `base_seed`, but the per-trial dummies use hardware entropy (selector seed $-1$) by design — required for a valid Monte Carlo FDR estimate. Re-runs reproduce the committed table to within MC noise ($\approx \pm 0.01$ at 200 trials).

---

## Output Files

Written to `simulation_results/` (already populated):

- `demo_trex_spca_01_mc_sim.txt` — human-readable table (methods × FDR/TPR/PEV × SNR column(s)).
- `demo_trex_spca_01_mc_sim.csv` — tidy long format: `method, metric, snr_db, value`.

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim
```

---

## Real Results (SNR = -10 dB, num_MC = 200, as committed)

| Method | FDR | TPR | PEV |
|---|---|---|---|
| OrdPCA | 0.9500 | 1.0000 | 0.2021 |
| OraclePCA | 0.0060 | 0.9940 | 0.1241 |
| TRexSPCA-EN-Act | 0.1229 | 0.9960 | 0.1216 |
| TRexSPCA-ENaug-Act | 0.1288 | 0.9960 | 0.1214 |
| TRexSPCA-EN-Thr | 0.1308 | 0.9970 | 0.1328 |
| TRexSPCA-ENaug-Thr | 0.1311 | 0.9970 | 0.1330 |

Reference points: the legacy CRAN R implementation reports FDR $0.100$ / TPR $0.999$ (Act) and FDR $0.100$ / TPR $1.000$ (Thr) at this same operating point. The R and Python counterparts of this demo (`R/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/`, `Python/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/`) reproduce the table above to within MC noise on their own RNG streams.

---

## Interpretation

- **OrdPCA**'s $\mathrm{FDR}\approx0.95$ and $\mathrm{TPR}=1.0$ are trivial artifacts of retaining all $p=100$ loadings (true support is $p_1=5$) — not a genuine comparison point, only a non-sparse reference.
- **OraclePCA** achieves near-perfect PC1 recovery (it is told the true support size); its PEV ($0.124$) is the useful reference value here — how much variance the true sparse signal explains, independent of any estimation procedure.
- All four **T-Rex SPCA** variants recover PC1's support almost perfectly (TPR $\approx 0.996$) with FDR close to the $\mathrm{tFDR}=0.10$ target — even though $-10$ dB sounds like a very weak signal, the strong first factor (std $5$) dominates PC1 recovery in this design on the covariance-PCA footing.
- The realized FDR ($\approx 0.12$–$0.13$) sits $\approx 0.02$–$0.03$ above the target and the legacy CRAN R reference ($0.100$). This is the long-known residual C++-vs-R gap investigated by the `validation/trex_spca/` programs (suspects: the internal CD-ridge λ₂ CV vs. R's `cv.glmnet`, and dummy-generation details); `validation_trex_spca_06_handrolled_comparison` confirms it is **not** a `TRexSPCA`-class artifact — the class and a hand-rolled per-PC pipeline replicating the legacy recipe agree with each other within MC noise.
- The four T-Rex variants are statistically indistinguishable from each other at this single SNR point; no solver (TENET vs. TENET_AUG) or mode (ActiveSet vs. Thresholded) conclusion should be drawn from it. The **Thresholded** variants' slightly higher PEV ($\approx 0.133$ vs. $\approx 0.122$) is the one systematic difference: they rebuild loadings from the full ordinary-PCA solution on the selected support, which preserves more explained variance than the active-set coefficients.
- If you extend this demo to a wider SNR grid, expect all methods to improve monotonically with SNR (the legacy reference reports T-Rex FDR falling to $\approx 0.06$ at $-7$ dB and $-5$ dB with TPR $= 1.0$) — treat exact values as hypotheses until run/committed.

---

**Last updated**: 2026-07-08
