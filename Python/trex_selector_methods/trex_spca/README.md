# T-Rex SPCA: Sparse Principal Component Analysis — Python Demos

## Overview

This folder contains Python example programs for **T-Rex SPCA**, a sparse principal component analysis method built on top of **T-Rex+GVS** (see [../trex_gvs/](../trex_gvs/README.md)). Instead of selecting sparse *variables* in a regression model, T-Rex SPCA selects sparse *loadings* for each principal component, applying the T-Rex+GVS elastic-net machinery per component to control the false discovery rate of the estimated loading support.

The single demo in this folder compares T-Rex SPCA against two PCA baselines on a synthetic sparse-factor-model dataset. It is the Python counterpart of the C++ suite in `cpp/trex_selector_methods/trex_spca/` and the R suite in `R/trex_selector_methods/trex_spca/`, using the `trex_selector_neo` `TRexSPCASelector`.

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

All methods are evaluated on the same $\boldsymbol{X}$ after a shared **center-only** preprocessing step (mean subtraction, no column scaling — a covariance-PCA footing), so that OrdPCA, OraclePCA, and T-Rex SPCA are compared on equal preprocessing. The column scales must **not** be normalized: in this factor model the amplitude signal lives in the column variances, and z-scoring the columns (correlation PCA) destroys it — measured effect at $-10$ dB: T-Rex SPCA FDR degrades from $\approx 0.13$ to $\approx 0.53$ (see the C++ `validation_trex_spca_06_handrolled_comparison`).

---

## T-Rex SPCA specific concepts

- **`ts.SPCAMode`** — how the sparse loading vector for each principal component is assembled from the underlying T-Rex+GVS selection:
  - `ActiveSet`: use only the T-Rex+GVS *selected* variables' loadings.
  - `Thresholded`: assemble loadings from the full ordinary-PCA solution, thresholded down to the T-Rex+GVS-selected support.
- **`ts.ENSolverType`** — the elastic-net solver driving each per-component T-Rex+GVS selection: `TENET` (Gram-based) or `TENET_AUG` (augmented-LASSO formulation); see [../trex_gvs/README.md](../trex_gvs/README.md) for the shared elastic-net background.
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
  ├── trex_spca_sim_utils.py
  └── demo_trex_spca_01_mc_sim/
      ├── demo_trex_spca_01_mc_sim.py
      ├── README.md
      └── simulation_results/
          ├── demo_trex_spca_01_mc_sim.txt
          └── demo_trex_spca_01_mc_sim.csv
```

The layout mirrors the C++ suite one-to-one. [trex_spca_sim_utils.py](trex_spca_sim_utils.py) is the suite-level shared module (Python counterpart of the C++ `trex_spca_sim_utils.hpp`): `generate_sparse_factor_model()` (numpy `default_rng`), `center_columns()`, `evaluate_spca()`, the parallel MC loops `run_mc_trials_pca()` / `run_mc_trials_oracle_pca()` / `run_mc_trials_trex_spca()` (multiprocessing `spawn` pool, default `min(6, cores)` workers), and `save_and_print_spca_mc_results()` (console + `.txt` table plus tidy CSV `method,metric,snr_db,value`).

---

## Example results

This demo has **real committed output** — from a single SNR point ($-10$ dB), `num_MC = 200`, produced by the center-only pipeline (run 2026-07-08):

| Method | FDR | TPR | PEV |
|---|---|---|---|
| OrdPCA | 0.9500 | 1.0000 | 0.2014 |
| OraclePCA | 0.0080 | 0.9920 | 0.1233 |
| TRexSPCA-EN-Act | 0.1317 | 0.9970 | 0.1202 |
| TRexSPCA-ENaug-Act | 0.1350 | 0.9970 | 0.1203 |
| TRexSPCA-EN-Thr | 0.1359 | 0.9980 | 0.1335 |
| TRexSPCA-ENaug-Thr | 0.1354 | 0.9980 | 0.1331 |

### Interpretation

- **OrdPCA**'s FDR $\approx 0.95$ is the expected trivial artifact of selecting all $p=100$ loadings (true support is only $p_1=5$) — not a meaningful comparison point.
- **OraclePCA** achieves near-perfect recovery (it is told the true support size); its PEV ($0.123$) is the reference for how much variance the true sparse signal explains, independent of any selection procedure.
- All four **T-Rex SPCA** variants recover PC1's support almost perfectly (TPR $\approx 0.997$–$0.998$) with FDR close to the $\mathrm{tFDR}=0.10$ target. The realized FDR ($\approx 0.13$–$0.14$) sits $\approx 0.03$ above the target and the legacy CRAN R reference ($0.100$ at this point) — the long-known residual gap between the C++ core and the CRAN R implementation, investigated by the `cpp/trex_selector_methods/validation/trex_spca/` programs; `validation_trex_spca_06_handrolled_comparison` confirms it is not a `TRexSPCA`-class artifact.
- The four T-Rex variants are statistically indistinguishable from each other at this one data point; no solver/mode conclusion should be drawn from a single SNR level. The **Thresholded** variants' slightly higher PEV ($\approx 0.133$ vs. $\approx 0.120$) is the one systematic difference.

> **History note (z-scoring episode, fixed 2026-07-08)**: an earlier version of this pipeline z-score-standardized $\boldsymbol{X}$ (and set the per-PC selector to `ScalingMode.ZSCORE`), which destroyed the factor amplitude signal and inflated the measured T-Rex FDR to $\approx 0.53$ at this same nominal SNR; the C++ `validation_trex_spca_06_handrolled_comparison` isolated the cause and the pipeline is back to center-only. `main()` currently exercises a single SNR value ($-10$ dB, explicitly commented `# TEMP fast-validation`) with `num_MC=200` — mirroring the C++ `main()`; the committed output matches this configuration. Trial-for-trial reproduction is intentionally impossible: each trial's **data** is drawn from a deterministic seed (`base_seed_offset + mc*1000`, with `base_seed_offset = base_seed + round(snr)`), but the per-trial dummies use hardware entropy (selector seed $-1$), as required for a valid Monte Carlo FDR estimate; re-runs reproduce the table to within MC noise ($\approx \pm 0.01$ at 200 trials). See the demo's own [README.md](demo_trex_spca_01_mc_sim/README.md) for detail, including the cross-language comparison against the committed C++ and R numbers.

---

## Running the demo

```bash
cd Python/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim
python demo_trex_spca_01_mc_sim.py                # default min(6, cores) workers
python demo_trex_spca_01_mc_sim.py 8              # 8 worker cores
```

Use the repo virtual environment (`.venv/bin/python`) if `trex_selector_neo` is not on the system interpreter. The demo inserts both its own directory and the parent suite directory onto `sys.path`, so the shared `trex_spca_sim_utils` module resolves from the nested location and is re-importable by the `spawn`-launched Monte Carlo workers — run it **as a file** (not via `python -` / piped stdin), matching the other Python suites. The optional first CLI argument sets the number of parallel worker cores.

Indices are **0-based** (Python convention); the R counterpart is 1-based.

---

## Notes for new users

- Start with `demo_trex_spca_01_mc_sim/` — it is currently the only demo in this folder.
- Read the PC1-only FDR/TPR rationale above before interpreting results — this is a deliberate evaluation-scope choice, not a bug.
- Cross-check numbers against the language counterparts: `cpp/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/` and `R/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/` run the identical simulation (on their own RNG streams — exact per-trial matches across languages are impossible: numpy PCG64 vs. C++ `std::mt19937`). The R cross-check probes live under `R/trex_selector_methods/validation/trex_spca/`.
- All three languages run the per-PC selector with the library default scaling (unit-L2, matching the legacy `lm_dummy.R` centering + L2 normalization); R's `trex_spca_control()` does not expose `scaling_mode`, but that is moot here.
- For deeper correctness checks of individual pipeline stages (λ₂ CV, solver equivalence, scaling comparison, step-by-step C++/R diffs), see the validation programs in `cpp/trex_selector_methods/validation/trex_spca/`.

---

## References

- See [../README.md](../README.md) for the category overview of all T-Rex selector variants.
- See [../trex_gvs/README.md](../trex_gvs/README.md) for the shared T-Rex+GVS elastic-net background (`ENSolverType`, `LambdaSelectionMethod`).

---

**Last updated**: 2026-07-08
