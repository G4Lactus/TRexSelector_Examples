# Demo 01: T-Rex SPCA Monte Carlo Simulation (Python)

## Purpose

Compare **T-Rex SPCA** (4 solver/mode combinations) against **ordinary PCA** and **oracle-thresholded PCA** on a synthetic sparse 3-factor model, evaluating PC1 loading-support recovery (FDR/TPR) and cumulative explained variance (PEV).

Python counterpart of `cpp/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim.cpp`. The R counterpart is `R/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/`. Support indices here are **0-based** (Python convention).

---

## Data Generation Parameters (`generate_sparse_factor_model`)

$$
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E}.
$$

- $\boldsymbol{Z}$ ($n \times M$): factor columns $\mathcal{N}(0, \sigma_m^2)$ with $\sigma_m \in \{5, 3, 1\}$ for $m=0,1,2$.
- $\boldsymbol{V}$ ($p \times M$): each column has exactly `p1` nonzero entries (value $0.9$), drawn without replacement from a shared candidate pool of size `overlap_pool_size` — factor supports can partially overlap.
- $\boldsymbol{E}$: i.i.d. Gaussian, scaled to hit a target SNR in dB.
- $n=50$, $p=100$, $p_1=5$, $M=3$, `overlap_pool_size=30`. RNG: numpy `default_rng` (PCG64).

All methods share one **center-only** $\boldsymbol{X}$ (see `center_columns()` in `../trex_spca_sim_utils.py`: mean subtraction, no column scaling). This puts OrdPCA, OraclePCA, and TRexSPCA on a common covariance-PCA footing, matching the legacy CRAN R reference. The column scales must **not** be normalized — the factor amplitude signal lives in the column variances, and z-scoring destroys it (see the History note below).

---

## Control Parameters

```
tFDR = 0.10
lambda_2 = -1        # < 0 triggers k-fold CV for ridge-penalty selection (0 = no ridge, > 0 = fixed)
num_MC = 200          # as set in main(), see caveat below
base_seed = 42        # seeds the per-trial DATA DGP deterministically:
                      # trial mc uses base_seed_offset + mc*1000,
                      # with base_seed_offset = base_seed + round(snr_db)
```

The T-Rex SPCA control is assembled as `ts.TRexSPCAControlParameter()` with `mode` (`ts.SPCAMode.ActiveSet` / `Thresholded`), `en_solver` (`ts.ENSolverType.TENET` / `TENET_AUG`), and the nested `gvs_ctrl` (`lambda2_method = ts.LambdaSelectionMethod.CV_1SE_CCD`, `lambda_2 = -1.0`; the per-PC selector scaling is left at the library default, unit-L2). The selector is constructed as `ts.TRexSPCASelector(np.asfortranarray(X), M, tFDR, ctrl, -1, False)`; `sel.select()` returns `None` and `res = sel.getResult()` exposes `res.V`, `res.Z`, and `res.active_sets` (0-based). The PCA baselines use the library `trex_selector_neo.ml_methods.PCA(center=False, normalize=False).fit(X, M)` → `PCAResult` with `.V` / `.Z` — `center=False` because $X$ is already centered (covariance PCA).

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

All four T-Rex SPCA variants run their per-component selector with the default `ScalingMode.L2` (matching the legacy `lm_dummy.R` column centering + L2 normalization), exactly as in C++ and R.

---

## Metrics

- **FDR / TPR** — evaluated strictly on **PC1's** loading support only; PCs 2–3 are excluded because ordinary PCA's orthogonality constraint mixes their supports across the true factors, making a per-component ground-truth comparison ambiguous beyond PC1.
- **PEV** — cumulative percentage of explained variance across all $M$ components, via QR decomposition of the estimated score matrix.

---

## ⚠️ Configuration notes

`main()` sets `num_MC = 200` and `snr_values = [-10.0]` — a **single SNR point**, explicitly commented `# TEMP fast-validation`, mirroring the C++ `main()`. Only this one SNR point has been run and committed; the file's doc-comment describes exactly this configuration.

> **History note (z-scoring episode, fixed 2026-07-08)**: an earlier version of this pipeline z-score-standardized $\boldsymbol{X}$ and set the per-PC selector to `ScalingMode.ZSCORE`. The z-scoring converts the covariance PCA into a correlation PCA and destroys the factor amplitude signal — the measured T-Rex FDR inflated to $\approx 0.52$ and OraclePCA TPR dropped to $\approx 0.86$ at this same nominal SNR. The C++ validation program `validation_trex_spca_06_handrolled_comparison` isolated the cause (the `TRexSPCA` class and a hand-rolled pipeline replicating the legacy R recipe agree within MC noise; only the preprocessing differed), and the pipeline is back to center-only.

Trial-for-trial reproduction is intentionally impossible: each trial's **data** is drawn from a deterministic seed (`base_seed_offset + mc*1000`), so the sequence of design matrices is fixed by `base_seed`, but the per-trial dummies use hardware entropy (selector seed $-1$) by design — required for a valid Monte Carlo FDR estimate. Re-runs reproduce the committed table to within MC noise ($\approx \pm 0.01$ at 200 trials).

---

## Output Files

Written to `simulation_results/` (already populated):

- `demo_trex_spca_01_mc_sim.txt` — human-readable table (methods × FDR/TPR/PEV × SNR column(s)).
- `demo_trex_spca_01_mc_sim.csv` — tidy long format: `method, metric, snr_db, value`.

---

## Running the Demo

```bash
cd Python/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim
python demo_trex_spca_01_mc_sim.py                # default min(6, cores) workers
python demo_trex_spca_01_mc_sim.py 8              # 8 worker cores
```

Use the repo virtual environment (`.venv/bin/python`) if `trex_selector_neo` is not on the system interpreter. The Monte Carlo trials run in parallel via a multiprocessing `spawn` pool, so the demo must be run **as a real `.py` file** (workers re-import the module) — not via `python -` / piped stdin. This matches the convention of the other Python suites. The optional first CLI argument sets the number of worker cores.

---

## Real Results (SNR = -10 dB, num_MC = 200, as committed, run 2026-07-08)

| Method | FDR | TPR | PEV |
|---|---|---|---|
| OrdPCA | 0.9500 | 1.0000 | 0.2014 |
| OraclePCA | 0.0080 | 0.9920 | 0.1233 |
| TRexSPCA-EN-Act | 0.1317 | 0.9970 | 0.1202 |
| TRexSPCA-ENaug-Act | 0.1350 | 0.9970 | 0.1203 |
| TRexSPCA-EN-Thr | 0.1359 | 0.9980 | 0.1335 |
| TRexSPCA-ENaug-Thr | 0.1354 | 0.9980 | 0.1331 |

Reference points: the legacy CRAN R implementation reports FDR $0.100$ / TPR $\approx 1.0$ at this same operating point. The committed C++ reference (same configuration, C++ RNG streams): OrdPCA $0.9500/1.0000/0.2021$; OraclePCA $0.0060/0.9940/0.1241$; T-Rex FDR $\approx 0.123$–$0.131$, TPR $\approx 0.996$–$0.997$; R: OraclePCA $0.0040/0.9960$, T-Rex FDR $\approx 0.129$–$0.141$, TPR $\approx 0.999$–$1.0$. Python agrees within MC noise — exact numbers cannot match because the data RNG streams differ across languages (numpy PCG64 vs. C++ `std::mt19937`).

---

## Interpretation

- **OrdPCA**'s $\mathrm{FDR}\approx0.95$ and $\mathrm{TPR}=1.0$ are trivial artifacts of retaining all $p=100$ loadings (true support is $p_1=5$) — not a genuine comparison point, only a non-sparse reference.
- **OraclePCA** achieves near-perfect PC1 recovery (it is told the true support size); its PEV ($0.123$) is the useful reference value here — how much variance the true sparse signal explains, independent of any estimation procedure.
- All four **T-Rex SPCA** variants recover PC1's support almost perfectly (TPR $\approx 0.997$–$0.998$) with FDR close to the $\mathrm{tFDR}=0.10$ target — even though $-10$ dB sounds like a very weak signal, the strong first factor (std $5$) dominates PC1 recovery on the covariance-PCA footing.
- The realized FDR ($\approx 0.13$–$0.14$) sits $\approx 0.03$ above the target and the legacy CRAN R reference ($0.100$). This is the long-known residual gap between the C++ core and the CRAN R implementation, investigated by the cross-check programs now in the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`) (suspects: the internal CD-ridge λ₂ CV vs. R's `cv.glmnet`, and dummy-generation details); `validation_trex_spca_06` confirms it is **not** a `TRexSPCASelector` artifact.
- The four T-Rex variants are statistically indistinguishable from each other at this single SNR point; no solver (TENET vs. TENET_AUG) or mode (ActiveSet vs. Thresholded) conclusion should be drawn from it. The **Thresholded** variants' slightly higher PEV ($\approx 0.133$ vs. $\approx 0.120$) is the one systematic difference: they rebuild loadings from the full ordinary-PCA solution on the selected support, which preserves more explained variance than the active-set coefficients.
- If you extend this demo to a wider SNR grid, expect all methods to improve monotonically with SNR (the legacy reference reports T-Rex FDR falling to $\approx 0.06$ at $-7$ dB and $-5$ dB with TPR $= 1.0$) — treat exact values as hypotheses until run/committed.

---

**Last updated**: 2026-07-08
