# T-Rex SPCA: Sparse PCA — R Demos and Cross-Check Probes

> **Status**: NEW — migrated to the **TRexSelectorNeo** R6 API. `demo_01` uses
> `TRexSPCASelector`; `demo_02` and the dummy-variance probes use
> `TRexGVSSelector`. No longer depends on the CRAN `TRexSelector` package. The
> two `glmnet` lambda2 probes are unchanged.

## Purpose

R programs for the **T-Rex Sparse PCA** method: a Monte Carlo evaluation
against ordinary PCA and an oracle baseline, an R reference-dump generator for
the C++/R head-to-head, and four probe scripts that cross-check CSV dumps
produced by the C++ validation demos.

---

## Demos

| File | Description |
|---|---|
| [demo_trex_spca_01.R](demo_trex_spca_01.R) | MC evaluation of T-Rex SPCA (via `TRexSPCASelector`) vs ordinary PCA, an oracle, and the CRAN `elasticnet::spca` baseline on a sparse M-factor model (n = 50, p = 100, M = 3 factors with overlapping support, SNR sweep −10 to +10 dB); reports TPR/FDR on PC1 support recovery and QR-adjusted cumulative explained variance (PEV) |
| [demo_trex_spca_02.R](demo_trex_spca_02.R) | R reference dump for the C++ head-to-head (via `TRexGVSSelector`): regenerates each dataset with the exact DGP of demo 01, dumps `X_<mc>.csv` at full precision plus `truth.csv`, `r_lambda2.csv`, `r_results.csv`, `r_pc1.csv` into `rdump/`, so the C++ pipeline can rerun on bit-identical inputs (uses `glmnet`) |

Console summaries of past runs are stored in
[simulation_results/](simulation_results/)
(`res_demo_trex_spca_01_x.txt`, `res_demo_trex_spca_02.txt`).

---

## Cross-check probes

These four scripts investigate the residual C++-vs-R FDR gap on the SPCA
benchmark. The first two **require the C++ validation demos to be run first**:
they read the CSV dumps (`lambda2_probe_X.csv`, `lambda2_probe_y.csv`) that
the C++ side writes into
`cpp/trex_selector_methods/validation/trex_spca/validation_results/` and
cross-check them with `glmnet`.

| File | Description |
|---|---|
| [lambda2_probe.R](lambda2_probe.R) | Loads the exact centered X and PC1 score y dumped by the C++ lambda2 probe; repeats `cv.glmnet` over several seeds and compares the mean LARS-unit ridge penalty `lambda_2_lars = lambda.1se * p / 2` against the C++ probe's `glmnet_lars` column |
| [lambda2_foldmatch.R](lambda2_foldmatch.R) | Decisive fold-matched head-to-head for `lambda.1se`: fixes a deterministic `foldid`, runs `cv.glmnet(alpha = 0)`, and dumps the full CV curve (lambda grid, cvm, cvsd) for point-by-point comparison against the C++ CD ridge CV on identical folds |
| [probe_R_dummy_variance.R](probe_R_dummy_variance.R) | Re-runs the R selector (`TRexGVSSelector`) on the shipped `rdump10/` X files under several global RNG seeds to test whether the −10 dB FDR value is systematic or dummy-draw noise (reads local CSVs; no fresh C++ run needed) |
| [r_dummy_variance_probe.R](r_dummy_variance_probe.R) | Quantifies the R selector's (`TRexGVSSelector`) own dummy-RNG variance on the identical dumped X in `rdump/` (reuses `X_*.csv`, `r_lambda2.csv`, `truth.csv`); compares the resulting FDR spread to the C++ band |

Data folders shipped with the repo: `rdump/`, `rdump10/`, and
`rdump_snr7_backup/` contain the per-trial `X_<mc>.csv` designs plus the R and
C++ result CSVs from previous head-to-head runs. Indices in these dump files
are **0-based** (converted to 1-based inside the R scripts).

Note: the probe headers refer to the C++ programs by older names
(`demo_trex_spca_02_lambda2_probe.cpp`, `demo_trex_spca_05_rdump_pipeline`,
`demo_trex_spca_06_lambda2_foldmatch.cpp`); these now live under
`validation/trex_spca/` as `validation_trex_spca_01_lambda2_probe.cpp`,
`validation_trex_spca_04_rdump_pipeline.cpp`, and
`validation_trex_spca_05_lambda2_foldmatch.cpp`.

---

## Running

```bash
Rscript R/trex_selector_methods/trex_spca/demo_trex_spca_01.R
Rscript R/trex_selector_methods/trex_spca/demo_trex_spca_02.R
# probes (run the corresponding cpp validation demo first where noted):
Rscript R/trex_selector_methods/trex_spca/lambda2_probe.R
```

Dependencies: `TRexSelectorNeo`, `elasticnet` (demo 01), `glmnet` (demo 02 and
the lambda2 probes), `parallel`.

---

## Counterparts

- C++ demo: [cpp/trex_selector_methods/trex_spca/](../../../cpp/trex_selector_methods/trex_spca/)
  — `demo_trex_spca_01_mc_sim`.
- C++ validation programs:
  [cpp/trex_selector_methods/validation/trex_spca/](../../../cpp/trex_selector_methods/validation/trex_spca/)
  — lambda2 probe, solver comparison, scaling comparison, rdump pipeline, and
  lambda2 foldmatch.

---

**Last updated**: 2026-07-08
