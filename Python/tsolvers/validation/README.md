# Validation: T-Algorithm Solver Correctness Checks (Python)

## Overview

Head-to-head numerical equivalence checks for individual solvers via the Python
bindings, comparing against reference implementations (R's CRAN `tlars` package,
or alternate formulations of the same estimator). Python mirror of
[cpp/tsolvers/validation/](../../../cpp/tsolvers/validation/README.md).

---

## Validations

| # | Name | Purpose |
|---|------|---------|
| 01 | [validation_ts_01_tenet_aug_comparison/](validation_ts_01_tenet_aug_comparison/README.md) | Three-way elastic-net equivalence: Gram-based `TENET_Solver` vs. augmented-LASSO `TENETAug_Solver` vs. manual `lasso_star` augmentation solved with `TLASSO_Solver`. Two scaling scenarios (Z-score, L2), complete paths; verdict: EQUIVALENT at machine precision. |
| 02 | [validation_ts_02_tlars_tlasso_rcompare/](validation_ts_02_tlars_tlasso_rcompare/README.md) | CRAN `tlars` vs. Python-binding `TLARS`/`TLASSO`/`TENET`/`TENETAug` path equivalence on SPCA sparse factor-model data. Reads the reference CSVs dumped by the R generator; mirrors the cpp comparator, supports `--dir` override, and gates with exit code 0/1. |

Validation 02 does **not** generate its own data: it consumes the reference CSVs
in `R/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/rdump_tlars/`,
produced by the CRAN-tlars generator script
`demo_ts_compare_tlars_tlasso.R` in the same R folder.

---

## Running

```bash
.venv/bin/python Python/tsolvers/validation/validation_ts_01_tenet_aug_comparison/validation_ts_01_tenet_aug_comparison.py
.venv/bin/python Python/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/validation_ts_02_tlars_tlasso_rcompare.py
```

---

**Last updated**: 2026-07-08
