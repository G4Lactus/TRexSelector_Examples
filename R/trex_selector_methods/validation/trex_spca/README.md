# Validation — T-Rex SPCA: R Cross-Check Programs

## Purpose

R-side counterparts of the C++ validation programs in
[cpp/trex_selector_methods/validation/trex_spca/](../../../../cpp/trex_selector_methods/validation/trex_spca/):
an R reference-dump generator for the C++/R head-to-head plus four probe
scripts that cross-check CSV dumps produced by the C++ validation demos.

These programs were relocated here from `R/trex_selector_methods/trex_spca/`
when that suite was restructured to mirror the cpp demo layout (2026-07-08);
the demo suite now contains only `demo_trex_spca_01_mc_sim/`.

---

## Programs

| File | Description |
|---|---|
| [demo_trex_spca_02.R](demo_trex_spca_02.R) | R reference dump for the C++ head-to-head (via `TRexGVSSelector`): regenerates each dataset with the exact DGP of the legacy R demo (X centered, **not** z-scored), dumps `X_<mc>.csv` at full precision plus `truth.csv`, `r_lambda2.csv`, `r_results.csv`, `r_pc1.csv` into `rdump/`, so the C++ pipeline can rerun on bit-identical inputs (uses `glmnet`) |
| [lambda2_probe.R](lambda2_probe.R) | Loads the exact centered X and PC1 score y dumped by the C++ lambda2 probe; repeats `cv.glmnet` over several seeds and compares the mean LARS-unit ridge penalty `lambda_2_lars = lambda.1se * p / 2` against the C++ probe's `glmnet_lars` column |
| [lambda2_foldmatch.R](lambda2_foldmatch.R) | Decisive fold-matched head-to-head for `lambda.1se`: fixes a deterministic `foldid`, runs `cv.glmnet(alpha = 0)`, and dumps the full CV curve (lambda grid, cvm, cvsd) for point-by-point comparison against the C++ CD ridge CV on identical folds |
| [probe_R_dummy_variance.R](probe_R_dummy_variance.R) | Re-runs the R selector (`TRexGVSSelector`) on the shipped `rdump10/` X files under several global RNG seeds to test whether the −10 dB FDR value is systematic or dummy-draw noise (reads local CSVs; no fresh C++ run needed) |
| [r_dummy_variance_probe.R](r_dummy_variance_probe.R) | Quantifies the R selector's (`TRexGVSSelector`) own dummy-RNG variance on the identical dumped X in `rdump/` (reuses `X_*.csv`, `r_lambda2.csv`, `truth.csv`); compares the resulting FDR spread to the C++ band |

`lambda2_probe.R` and `lambda2_foldmatch.R` **require the C++ validation demos
to be run first**: they read the CSV dumps (`lambda2_probe_X.csv`,
`lambda2_probe_y.csv`) that the C++ side writes into
`cpp/trex_selector_methods/validation/trex_spca/validation_results/`.

The recorded console output of a past `demo_trex_spca_02.R` run is kept in
[res_demo_trex_spca_02.txt](res_demo_trex_spca_02.txt).

---

## Data folders

`rdump/`, `rdump10/`, and `rdump_snr7_backup/` contain the per-trial
`X_<mc>.csv` designs plus the R and C++ result CSVs from previous head-to-head
runs. Indices in these dump files are **0-based** (converted to 1-based inside
the R scripts).

The C++ program `validation_trex_spca_04_rdump_pipeline` reads/writes `rdump/`
at this location by default
(`<repo>/R/trex_selector_methods/validation/trex_spca/rdump`); override with
`--dir <path>`.

Note: the probe headers refer to the C++ programs by older names
(`demo_trex_spca_02_lambda2_probe.cpp`, `demo_trex_spca_05_rdump_pipeline`,
`demo_trex_spca_06_lambda2_foldmatch.cpp`); these now live under
`validation/trex_spca/` as `validation_trex_spca_01_lambda2_probe.cpp`,
`validation_trex_spca_04_rdump_pipeline.cpp`, and
`validation_trex_spca_05_lambda2_foldmatch.cpp`.

---

## Running

```bash
Rscript R/trex_selector_methods/validation/trex_spca/demo_trex_spca_02.R
# probes (run the corresponding cpp validation demo first where noted):
Rscript R/trex_selector_methods/validation/trex_spca/lambda2_probe.R
```

Dependencies: `TRexSelectorNeo`, `glmnet` (demo 02 and the lambda2 probes).

---

## Counterparts

- Demo suite: [R/trex_selector_methods/trex_spca/](../../trex_spca/README.md)
  — `demo_trex_spca_01_mc_sim`.
- C++ validation programs:
  [cpp/trex_selector_methods/validation/trex_spca/](../../../../cpp/trex_selector_methods/validation/trex_spca/)
  — lambda2 probe, solver comparison, scaling comparison, rdump pipeline, and
  lambda2 foldmatch.

---

**Last updated**: 2026-07-08
