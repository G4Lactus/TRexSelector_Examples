# Validation: T-Algorithm Solver Correctness Checks (R)

## Overview

Head-to-head numerical correctness material for the T-Algorithm solvers on the
R side. This folder plays two roles: it ports the solver equivalence check to
the R6 bindings, and it hosts the **CRAN-tlars reference generator** whose
dumps drive the cross-language comparators in cpp and Python. Counterpart of
[cpp/tsolvers/validation/](../../../cpp/tsolvers/validation/README.md).

---

## Validations

| # | Name | Purpose |
|---|------|---------|
| 01 | [validation_ts_01_tenet_aug_comparison/](validation_ts_01_tenet_aug_comparison/README.md) | R port of the cpp program: Gram-based `TENET_Solver` vs. augmented-LASSO `TENETAug_Solver` vs. manual `lasso_star` augmentation via `TLASSO_Solver`, using the R6 solvers. The R binding has no `get_beta_path()`, so the X-block path matrix is assembled by looping `get_beta(step)` over the steps. |
| 02 | [validation_ts_02_tlars_tlasso_rcompare/](validation_ts_02_tlars_tlasso_rcompare/README.md) | CRAN-tlars reference **generator** (`demo_ts_compare_tlars_tlasso.R`) plus its `rdump_tlars/` CSV dumps. Consumed by BOTH the cpp comparator `validation_ts_02_tlars_tlasso_rcompare` and the Python comparator of the same name. |

---

## On the CRAN dependency in validation 02

The generator script deliberately uses the CRAN `tlars` package (and CRAN
`TRexSelector` internals for the GVS dummy-from-clustering reference) as an
**independent external ground truth** for the forward-selection
path-equivalence test of the native solvers. It is not legacy code awaiting
migration — replacing it with the bindings it validates would defeat the
purpose of the check.

---

## Running

```bash
Rscript R/tsolvers/validation/validation_ts_01_tenet_aug_comparison/validation_ts_01_tenet_aug_comparison.R
Rscript R/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/demo_ts_compare_tlars_tlasso.R   # requires CRAN tlars + TRexSelector
```

---

**Last updated**: 2026-07-08
