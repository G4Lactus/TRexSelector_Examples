# Validation: T-Algorithm Solver Correctness Checks

## Overview

Head-to-head numerical equivalence checks for individual solvers, comparing against reference implementations (R's `tlars` package, or alternate C++ formulations of the same estimator).

---

## Validations

| # | Name | Purpose |
|---|------|---------|
| 01 | [validation_ts_01_tenet_aug_comparison/](validation_ts_01_tenet_aug_comparison/README.md) | Gram-based `TENET_Solver` vs. augmented-LASSO `TENETAug_Solver` equivalence |
| 02 | [validation_ts_02_tlars_tlasso_rcompare/](validation_ts_02_tlars_tlasso_rcompare/README.md) | CRAN `tlars` vs. C++ `TLARS`/`TLASSO` path equivalence on SPCA data |

---

## Output Convention

Validation demos route their output (if any) to `validation_results/` folders instead of `simulation_results/`, distinguishing correctness checks from empirical Monte Carlo studies.

---

## Running

```bash
./build/debug/bin/tsolvers/validation/validation_ts_01_tenet_aug_comparison/validation_ts_01_tenet_aug_comparison
./build/debug/bin/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/validation_ts_02_tlars_tlasso_rcompare
```

---

**Last updated**: 2026-07-01
