# Validation: HAC Clustering Correctness Checks

## Overview

This folder contains validation programs for the hierarchical agglomerative clustering (HAC) components in `ml_methods`.

These checks are heavier than standard unit tests: they compare full C++ clustering workflows against R reference outputs on the same dumped design matrices and derived objects.

| # | Name | Purpose |
|---|------|---------|
| 01 | [validation_mlm_hac_01_rcompare/](validation_mlm_hac_01_rcompare/README.md) | Verifies that C++ single-linkage HAC reproduces R `hclust(method="single")` on identical datasets |
| 02 | [validation_mlm_hac_02_gvs_dummy_rcompare/](validation_mlm_hac_02_gvs_dummy_rcompare/README.md) | Verifies that C++ GVS dummy generation uses the same dendrogram cut and the same per-cluster generating covariance as R `TRexSelector:::add_dummies_GVS()` |

---

## Running

```bash
./build/debug/bin/ml_methods/validation/hac_clustering/validation_mlm_hac_01_rcompare/validation_mlm_hac_01_rcompare
./build/debug/bin/ml_methods/validation/hac_clustering/validation_mlm_hac_02_gvs_dummy_rcompare/validation_mlm_hac_02_gvs_dummy_rcompare
```

---

## Notes

Both validations expect reference CSV dumps produced by the corresponding R-side workflow.

---

**Last updated**: 2026-07-01
