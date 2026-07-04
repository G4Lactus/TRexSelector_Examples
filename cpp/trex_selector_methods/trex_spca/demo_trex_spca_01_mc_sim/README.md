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

All methods share one z-score-standardized $\boldsymbol{X}$ (see `standardize_columns()` — note: the current implementation centers columns and returns early, so variance normalization is effectively dead code / a left-in-place A/B test; all methods still see the same, centered-only, $\boldsymbol{X}$).

---

## Control Parameters

```
tFDR = 0.10
lambda_2 = -1        # triggers k-fold CV for ridge-penalty selection
num_MC = 200          # as set in main(), see caveat below
base_seed = 42
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

All four T-Rex SPCA variants run their per-component selector with `ScalingMode::ZSCORE`.

---

## Metrics

- **FDR / TPR** — evaluated strictly on **PC1's** loading support only; PCs 2–3 are excluded because ordinary PCA's orthogonality constraint mixes their supports across the true factors, making a per-component ground-truth comparison ambiguous beyond PC1.
- **PEV** — cumulative percentage of explained variance across all $M$ components, via QR decomposition of the estimated score matrix.

---

## ⚠️ Configuration / documentation mismatch

Three different numbers appear for this one demo, and they do not agree:

1. The file's own top-of-file doc-comment describes **"SNR sweep, `snr_db` in $\{-10,-5,0,5,10\}$ dB, num_MC=100"**.
2. The actual `main()` function sets `cfg.num_MC = 200` but `snr_values = {-10.0}` — a **single SNR point**, explicitly commented `// TEMP fast-validation`.
3. The **committed output file** (`simulation_results/demo_trex_spca_01_mc_sim.txt`) reports results **"averaged over 80 MC trials"**.

None of 100, 200, or 80 agree, and only one of the five documented SNR points has ever actually been run and committed. Treat the committed numbers below as illustrative of the $-10$ dB point only, generated under whatever configuration was active at the time — not as a reproducible guarantee if you rebuild and rerun today's `main()`.

---

## Output Files

Written to `simulation_results/` (already populated):

- `demo_trex_spca_01_mc_sim.txt` — human-readable table (methods × FDR/TPR/PEV × SNR column(s)).
- `demo_trex_spca_01_mc_sim.csv` — tidy long format: `method, metric, snr_db, value`.

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_spca_01_mc_sim
```

---

## Real Results (SNR = -10 dB, as committed)

| Method | FDR | TPR | PEV |
|---|---|---|---|
| OrdPCA | 0.9500 | 1.0000 | 0.2024 |
| OraclePCA | 0.0000 | 1.0000 | 0.1254 |
| TRexSPCA-EN-Act | 0.1188 | 0.9975 | 0.1204 |
| TRexSPCA-ENaug-Act | 0.1337 | 0.9975 | 0.1205 |
| TRexSPCA-EN-Thr | 0.1099 | 1.0000 | 0.1252 |
| TRexSPCA-ENaug-Thr | 0.1107 | 1.0000 | 0.1253 |

---

## Interpretation

- **OrdPCA**'s $\mathrm{FDR}\approx0.95$ and $\mathrm{TPR}=1.0$ are trivial artifacts of retaining all $p=100$ loadings (true support is $p_1=5$) — not a genuine comparison point, only a non-sparse reference.
- **OraclePCA** achieves perfect PC1 recovery by construction; its PEV ($0.1254$) is the useful reference value here — how much variance the true sparse signal explains, independent of any estimation procedure.
- All four **T-Rex SPCA** variants recover PC1's support almost perfectly (TPR $\approx 0.998$–$1.0$) with FDR close to the $\mathrm{tFDR}=0.10$ target, despite the nominally low $-10$ dB SNR — the strong first factor (std $5$) appears to dominate PC1 recovery in this design even when total signal-to-noise (across the whole matrix) is nominally weak.
- **Thresholded** mode is marginally ahead of **ActiveSet** mode here (slightly higher TPR, FDR slightly closer to target), and **TENET** is marginally ahead of **TENET_AUG** within each mode — both are single-point differences from one SNR level and MC configuration, not general conclusions; re-running additional SNR points (as the file's own doc-comment intends) would be needed to confirm these patterns hold more broadly.
- If you extend this demo to the full 5-point SNR grid described in the file's doc-comment, expect all four T-Rex SPCA variants' TPR to fall and FDR to become noisier as SNR decreases toward $-10$ dB and below, with OraclePCA's gap to the T-Rex variants widening — but this has not yet been run/committed, so treat it as a hypothesis to check, not a stated result.

---

**Last updated**: 2026-07-04
