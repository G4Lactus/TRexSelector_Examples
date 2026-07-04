# Validation: TENET vs. TENETAug Equivalence

## Purpose

Equivalence check between two elastic-net formulations:
- **`TENET_Solver`**: Gram-matrix-based elastic net
- **`TENETAug_Solver`**: augmented-LASSO formulation of elastic net (data augmentation trick)

Tested under two preprocessing scenarios (external vs. internal normalization) to confirm both solvers produce numerically identical coefficient paths.

---

## Running

```bash
./build/debug/bin/tsolvers/validation/validation_ts_01_tenet_aug_comparison/validation_ts_01_tenet_aug_comparison
```

---

## What to Look For

- Coefficient path agreement (within numerical tolerance) between `TENET_Solver` and `TENETAug_Solver`
- Consistency across both preprocessing scenarios

---

**Last updated**: 2026-07-01
