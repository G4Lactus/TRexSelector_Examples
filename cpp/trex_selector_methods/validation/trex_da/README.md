# Validation: DA-TRex Correctness Checks

## Overview

This folder contains validation programs for the Dependency-Aware T-Rex (DA-TRex) selector in `trex_selector_methods`.

These checks are heavier than standard unit tests: they compare full C++ selector/clustering workflows against an R reference implementation on the same fixed dataset.

| # | Name | Purpose |
|---|------|---------|
| 01 | [validation_trex_da_01_bt_dendro_diag/](validation_trex_da_01_bt_dendro_diag/README.md) | Verifies that the C++ DA-BT hierarchical clustering pipeline (dendrogram heights, per-$\rho$ group memberships, and final selection) reproduces an equivalent R implementation on one identical fixed dataset |

---

## Running

```bash
./build/debug/bin/trex_selector_methods/validation/trex_da/validation_trex_da_01_bt_dendro_diag/validation_trex_da_01_bt_dendro_diag
```

---

## Notes

Validation 01 is a single-dataset correctness diagnostic, not a Monte Carlo study — see its own README for the companion R script and PASS/FAIL interpretation.

---

**Last updated**: 2026-07-04
