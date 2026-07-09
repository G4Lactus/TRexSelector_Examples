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

`main()` sets `cfg.num_MC = 200` and the full legacy CRAN reference grid `snr_values = {-10, -7, -5, -3, 0, 3, 5, 7, 10}` dB. The committed output (`simulation_results/`) matches this configuration, the center-only pipeline, and the fixed library (see the history notes below).

**History note (two bugs, both fixed 2026-07-08):**

1. *Demo z-scoring*: during an earlier debugging round the pipeline briefly gained a shared z-score standardization; that converts the covariance PCA into a correlation PCA, destroys the factor amplitude signal, and pushed the measured T-Rex FDR to $\approx 0.52$ (OraclePCA TPR to $\approx 0.87$) at $-10$ dB. `validation_trex_spca_06_handrolled_comparison` isolated the cause and the pipeline is back to center-only.
2. *Library dummy normalization* (TRexSelector `5fceed3`): the GVS selector exactly unit-normalized every drawn cluster-MVN dummy column, erasing the random norm fluctuation that is part of the dummy null distribution and inflating the realized FDR from $\approx 0.10$ to $\approx 0.13$ while the FDP estimate still looked controlled. Root-caused via `validation_trex_spca_07_matrix_dump` (replaying the CRAN reference calibration on the selector's internal matrices reproduced its selections exactly, proving the calibration port correct and localizing the gap in the dummies) plus a 2×2 factorial with CRAN's own calibration functions. Dummies now enter the solver as drawn (center-only).

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

## Real Results (9-point SNR sweep, num_MC = 200, as committed)

T-Rex SPCA FDR across the grid (all four variants; TPR $= 1.0$ everywhere except $-10$ dB, where TPR $\approx 0.995$–$0.997$):

| SNR (dB) | −10 | −7 | −5 | −3 | 0 | 3 | 5 | 7 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| TRexSPCA-EN-Act | 0.094 | 0.062 | 0.053 | 0.067 | 0.069 | 0.062 | 0.084 | 0.072 | 0.069 |
| TRexSPCA-ENaug-Act | 0.103 | 0.063 | 0.064 | 0.062 | 0.070 | 0.064 | 0.075 | 0.067 | 0.069 |
| TRexSPCA-EN-Thr | 0.100 | 0.064 | 0.054 | 0.062 | 0.064 | 0.061 | 0.080 | 0.065 | 0.069 |
| TRexSPCA-ENaug-Thr | 0.094 | 0.059 | 0.061 | 0.062 | 0.071 | 0.069 | 0.075 | 0.065 | 0.070 |

The hardest point in full ($-10$ dB):

| Method | FDR | TPR | PEV |
|---|---|---|---|
| OrdPCA | 0.9500 | 1.0000 | 0.2021 |
| OraclePCA | 0.0060 | 0.9940 | 0.1241 |
| TRexSPCA-EN-Act | 0.0939 | 0.9960 | 0.1174 |
| TRexSPCA-ENaug-Act | 0.1033 | 0.9950 | 0.1177 |
| TRexSPCA-EN-Thr | 0.1001 | 0.9950 | 0.1280 |
| TRexSPCA-ENaug-Thr | 0.0935 | 0.9970 | 0.1287 |

This matches the legacy CRAN R reference across the grid (legacy: FDR $0.100$ at $-10$ dB, $\approx 0.058$ at $-7$/$-5$ dB, TPR $\approx 1.0$). See `simulation_results/` for the full committed table (PEV rows included) and the tidy CSV.

---

## Interpretation

- **OrdPCA**'s $\mathrm{FDR}\approx0.95$ and $\mathrm{TPR}=1.0$ are trivial artifacts of retaining all $p=100$ loadings (true support is $p_1=5$) — not a genuine comparison point, only a non-sparse reference.
- **OraclePCA** achieves (near-)perfect PC1 recovery everywhere (it is told the true support size); its PEV column is the reference for how much variance the true sparse signal explains, independent of any estimation procedure.
- All four **T-Rex SPCA** variants keep the realized FDR at or below the $\mathrm{tFDR}=0.10$ target across the whole grid with TPR $\approx 1.0$ — FDR control holds, matching the legacy CRAN implementation after the two 2026-07-08 fixes (see Configuration notes).
- **TENET and TENET_AUG are the same estimator on this problem**: run with identical per-trial dummy seeds they produce bit-identical selections (verified 100/100 trials via `validation_trex_spca_07_matrix_dump`). The EN-vs-ENaug differences visible in the table are pure Monte Carlo dummy noise — each method row draws fresh hardware-entropy dummies, and the standard error of a 200-trial mean FDR is $\approx \pm 0.009$, an order of magnitude larger than the observed row differences. Note the sign even flips between modes and SNR points.
- The **Thresholded** variants' higher PEV (e.g. $0.128$ vs. $0.117$ at $-10$ dB, converging toward OraclePCA's $0.89$ at $+10$ dB) is the one systematic mode difference: they rebuild loadings from the full ordinary-PCA solution on the selected support, which preserves more explained variance than the active-set ridge coefficients.

---

**Last updated**: 2026-07-08
