# Validation: TENET vs. TENETAug Equivalence (Python)

## Purpose

Equivalence check across three elastic-net computations on the same data,
via the Python bindings (mirrors the cpp validation of the same name):

- **`TENET_Solver`**: Gram-matrix-based elastic net
- **`TENETAug_Solver`**: augmented-LASSO formulation of elastic net (data augmentation trick)
- **manual `TLASSO_Solver` on hand-augmented matrices**: a literal replica of the R `lasso_star` reference (plain LASSO on the augmented `X_aug`/`y_aug` system, then de-normalised by `1/d2`)

All three are pre-scaled *externally* before the solvers run (each solver is
constructed with `normalize=False`, `intercept=False`), and the check is
repeated under two scaling scenarios:

1. **Z-score** scaling (centre + Bessel-corrected sample standard deviation), and
2. **L2-norm** scaling (centre + column L2 norm).

The goal is to confirm the Gram-based and augmented-LASSO elastic-net paths
agree with each other and with the R `lasso_star` reference.

---

## Configuration

- `n = 90`, `p = 150`, `L = 450` dummies, `lambda2 = 100.01`
- True support (0-based): `{10, 50, 85}` with coefficients `{2.5, -1.8, 3.2}`
- `snr = 1`, seed `42`
- Complete paths (`early_stop=False`; the nominal `T_stop` is unused)

---

## Running

```bash
.venv/bin/python Python/tsolvers/validation/validation_ts_01_tenet_aug_comparison/validation_ts_01_tenet_aug_comparison.py
```

---

## What to Look For

- Pairwise L2 agreement (within numerical tolerance) among the three paths: `TENET`, `TENETAug`, and the manual `lasso_star` replica
- Per-step RMSE to the true beta (de-normalised back to the original data scale)
- Consistency across both the Z-score and L2-norm scaling scenarios

**Verified result**: all three formulations are reported EQUIVALENT in both
scaling scenarios, with an overall worst pairwise max L2 of ~1e-12 (machine
precision).

---

**Last updated**: 2026-07-08
