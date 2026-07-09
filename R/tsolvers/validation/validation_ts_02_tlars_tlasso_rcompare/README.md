# Validation: CRAN `tlars` Reference Generator

## Purpose

Reference **generator** for the forward-selection path-equivalence test of the
native T-Algorithm solvers against an **independent external ground truth**:
the CRAN `tlars` package (plus CRAN `TRexSelector` internals for the GVS
dummy-from-clustering reference). This is deliberate — the script is not
legacy code awaiting migration; it exists precisely because it does *not* use
the bindings it validates.

The data-generating process is the same sparse factor-model DGP as the SPCA
demos: strongly correlated columns (factor structure + an overlap pool) that
actually trigger LARS ties and collinearity/LASSO sign-change drops.

Solver pairings covered by the dumps:

- **LARS** — CRAN `tlars` `type = "lar"` vs native `TLARS_Solver`
- **LASSO** — CRAN `tlars` `type = "lasso"` vs native `TLASSO_Solver`
- **EN** — CRAN augmented lasso (`lasso_star` recipe) with LASSO inner vs
  native `TENET_Solver` (Gram EN) and `TENETAug_Solver` (default inner), and
  with pure-LARS inner vs `TENETAug_Solver(use_lars_inner)`

Key API mapping: CRAN `tlars` takes one augmented matrix `Z = cbind(X, D)`
plus `num_dummies = L`, while the native solvers take X and D separately. To
remove preprocessing confounds, Z is pre-standardized (center + unit L2 per
column) and y centered once, then dumped verbatim so both sides see
bit-identical inputs, with the complete path computed (`early_stop = FALSE`).

---

## `rdump_tlars/` reference CSVs

The script writes its dumps into `rdump_tlars/` next to the script; the
comparators read them back:

| File(s) | Content |
|---|---|
| `meta.csv` | Single row: `n, p, L, num, lambda2, snr_db` |
| `Xn_<i>.csv`, `Dn_<i>.csv`, `y_<i>.csv` | Pre-standardized design (n x p), dummies (n x L), centered response, per trial |
| `r_lar_beta_<i>.csv` / `r_lar_act_<i>.csv` | LARS reference coefficient path + action sequence |
| `r_lasso_beta_<i>.csv` / `r_lasso_act_<i>.csv` | LASSO reference path + actions |
| `r_en_beta_<i>.csv` / `r_en_act_<i>.csv` | EN (augmented lasso, LASSO inner) reference path + actions |
| `r_enlar_beta_<i>.csv` / `r_enlar_act_<i>.csv` | EN (augmented lasso, pure-LARS inner) reference path + actions |
| `r_clust_height_<i>.csv` / `r_clust_labels_<i>.csv` | Single-linkage correlation-distance clustering reference (dendrogram heights + per-K cutree labels) |
| `gvs_corrmax.csv`, `r_gvs_labels_<i>_<cm>.csv`, `r_gvs_sigma_<i>_<cm>.csv` | GVS dummy-from-clustering reference per `corr_max` cut point: cluster labels + block-diagonal generating covariance from `TRexSelector:::add_dummies_GVS` |

---

## Running

Requires the CRAN packages `tlars` and `TRexSelector`.

```bash
Rscript R/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/demo_ts_compare_tlars_tlasso.R
```

---

## Consumers

The dumps are read by **two** comparators, which gate the native solvers
against the CRAN reference:

- C++: [cpp/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/](../../../../cpp/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/README.md)
- Python: [Python/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/](../../../../Python/tsolvers/validation/validation_ts_02_tlars_tlasso_rcompare/README.md)

---

**Last updated**: 2026-07-08
