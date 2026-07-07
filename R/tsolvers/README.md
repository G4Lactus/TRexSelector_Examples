# T-Solvers — R Cross-Validation Demo

## Purpose

Cross-validation of the C++ T-Algorithm solvers against an **independent
reference implementation**: the CRAN `tlars` package. This is deliberate — the
script uses CRAN `tlars` (and CRAN TRexSelector internals for the GVS
dummy-from-clustering reference) as an external ground truth for a
forward-selection **path-equivalence test** of the C++ solvers; it is not
legacy code awaiting migration.

---

## Demo

| File | Description |
|---|---|
| [demo_ts_compare_tlars_tlasso.R](demo_ts_compare_tlars_tlasso.R) | Reference generator for the TLARS / TLASSO / TENET path-equivalence test on the sparse factor-model DGP (strongly correlated columns that actually trigger LARS ties and LASSO sign-change drops) |

Solver pairings exercised:

- **LARS** — CRAN `tlars` `type = "lar"` vs C++ `TLARS_Solver`;
- **LASSO** — CRAN `tlars` `type = "lasso"` vs C++ `TLASSO_Solver`;
- **EN** — CRAN augmented-lasso ("lasso_star" recipe) vs C++ `TENET_Solver`
  (Gram-based EN) and `TENETAug_Solver` (augmented-lasso EN).

Key API mapping: CRAN `tlars` takes one augmented matrix
`Z = cbind(X, D)` plus `num_dummies = L`, while the C++ solvers take X and D
separately; R column `j` (1-based) corresponds to the C++ combined 0-based
index `j - 1`. To remove preprocessing confounds, Z is pre-standardized
(center + unit L2 per column) and y centered once, then dumped verbatim so
both sides see bit-identical inputs, with the complete path computed
(`early_stop = FALSE`).

---

## `rdump_tlars/` reference CSVs

The script writes its dumps into `rdump_tlars/` next to the script; the C++
comparison program reads them back:

| File(s) | Content |
|---|---|
| `meta.csv` | Single row: `n, p, L, num, lambda2, snr_db` |
| `Xn_<i>.csv`, `Dn_<i>.csv`, `y_<i>.csv` | Pre-standardized design (n x p), dummies (n x L), centered response, per trial |
| `r_lar_beta_<i>.csv` / `r_lar_act_<i>.csv` | LARS reference coefficient path + action sequence |
| `r_lasso_beta_<i>.csv` / `r_lasso_act_<i>.csv` | LASSO reference path + actions |
| `r_en_beta_<i>.csv` / `r_en_act_<i>.csv` | EN (augmented lasso) reference path + actions |

---

## Running

Requires the CRAN packages `tlars` and `TRexSelector`.

```bash
Rscript R/tsolvers/demo_ts_compare_tlars_tlasso.R
```

---

## Counterparts and coverage

- C++ comparison program:
  [cpp/tsolvers/validation/](../../cpp/tsolvers/validation/) —
  `validation_ts_02_tlars_tlasso_rcompare` (the script header still cites the
  older name `demo_ts_14_tlars_tlasso_rcompare.cpp`).
- The C++ tree ships **12 per-solver demos**
  (`demo_ts_01_tlars` … `demo_ts_12_tafs`, see
  [cpp/tsolvers/](../../cpp/tsolvers/)) that are not yet mirrored in R; this
  folder currently contains only the cross-validation script above.

---

**Last updated**: 2026-07-06
