# Validation: TENET vs. TENETAug Equivalence (R)

## Purpose

Equivalence check across three elastic-net computations on the same data,
via the R6 bindings (R port of the cpp validation of the same name):

- **`TENET_Solver`**: Gram-matrix-based elastic net
- **`TENETAug_Solver`**: augmented-LASSO formulation of elastic net (data augmentation trick)
- **manual `TLASSO_Solver` on hand-augmented matrices**: a literal replica of the R `lasso_star` reference (plain LASSO on the augmented `X_aug`/`y_aug` system, then de-normalised by `1/d2`)

All three are pre-scaled *externally* before the solvers run (each solver is
constructed with `normalize = FALSE`, `intercept = FALSE`), and the check is
repeated under two scaling scenarios:

1. **Z-score** scaling (centre + Bessel-corrected sample standard deviation), and
2. **L2-norm** scaling (centre + column L2 norm).

The goal is to confirm the Gram-based and augmented-LASSO elastic-net paths
agree with each other and with the R `lasso_star` reference.

Binding note: the R binding has no `get_beta_path()` accessor, so the p x
(steps+1) X-block path matrix is assembled by looping `get_beta(step)` over
steps `0..get_num_steps()` — the X block is exactly what the cpp comparison
uses.

---

## Configuration

- `n = 90`, `p = 150`, `L = 450` dummies, `lambda2 = 100.01`
- True support (1-based): `{11, 51, 86}` with coefficients `{2.5, -1.8, 3.2}`
  (the cpp program's 0-based `{10, 50, 85}`)
- `snr = 1`, seed `42`
- Complete paths (`early_stop = FALSE`; the nominal `T_stop` is unused)

---

## Running

```bash
Rscript R/tsolvers/validation/validation_ts_01_tenet_aug_comparison/validation_ts_01_tenet_aug_comparison.R
```

---

## What to Look For

- Pairwise L2 agreement (within numerical tolerance) among the three paths: `TENET`, `TENETAug`, and the manual `lasso_star` replica
- Per-step RMSE to the true beta (de-normalised back to the original data scale)
- Consistency across both the Z-score and L2-norm scaling scenarios

**Verified result**: all three formulations are reported EQUIVALENT in both
scaling scenarios, with an overall worst pairwise max L2 of ~8e-13 (machine
precision).

---

**Last updated**: 2026-07-08
