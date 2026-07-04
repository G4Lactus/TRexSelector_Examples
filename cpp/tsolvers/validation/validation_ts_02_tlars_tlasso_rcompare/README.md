# Validation: T-LARS/T-LASSO vs. R `tlars` Path Equivalence

## Purpose

Final forward-selection path-equivalence test comparing the CRAN `tlars` package against the C++ `TLARS`/`TLASSO` solvers, using the same sparse factor-model data as the SPCA demos (`trex_spca`).

This is the authoritative cross-language correctness check for the two most commonly used solvers in the T-Rex selector.

---

## Running

```bash
./build/debug/bin/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/validation_ts_02_tlars_tlasso_rcompare
```

---

## What to Look For

- Exact (or near machine-precision) agreement between C++ and R selection paths and coefficients
- Any divergence points (e.g., tie-breaking differences) should be flagged and documented

---

**Last updated**: 2026-07-01
