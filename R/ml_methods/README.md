# ML Methods — R Demos

Data-preprocessing and numerical building blocks used by the T-Rex solver and
selector layers, demonstrated through the **TRexSelectorNeo** R API (plus one
glmnet-based cross-check tool).

---

## Contents

| Folder | Status | What it covers |
|---|---|---|
| [hac_clustering/](hac_clustering/README.md) | NEW (TRexSelectorNeo) | Agglomerative hierarchical clustering via `agglomerative_cluster()`, in-memory and mmap-backed, with ARI evaluation against known ground truth |
| [standardization/](standardization/README.md) | NEW (TRexSelectorNeo) | Column standardization via the `ZScoreScaler` and `LpNormScaler` R6 classes, in-memory and mmap-backed |
| [pca/](pca/README.md) | NEW (TRexSelectorNeo) | Principal component analysis via the `PCA` R6 class: fit, restore round-trip, out-of-sample transform |
| [svd/](svd/README.md) | NEW (TRexSelectorNeo) | Truncated singular value decomposition via `SVDSolver` (direct and Gram paths, orthogonality checks) |
| [ridge_regression/](ridge_regression/README.md) | NEW (TRexSelectorNeo) | L2-penalized least squares via `RidgeSolver` (primal and dual paths, lambda sweep) |
| [model_selection/](model_selection/README.md) | NEW (TRexSelectorNeo) + cross-check (glmnet) | K-fold cross-validated ridge (`RidgeCV`) and elastic net (`ElasticNetCV`) with lambda.min / lambda.1se, plus the glmnet reference generator for the C++ coordinate-descent elastic net |

---

## Coverage relative to the C++ tree

The R tree now mirrors every [cpp/ml_methods/](../../cpp/ml_methods/) demo
area. Remaining difference:

- The R `ElasticNetCV` binding exposes cross-validation but no standalone
  elastic-net *path* class, so the R elastic-net demo covers the CV scenarios
  only; the full path fit (cpp `ElasticNet.fit`/`fit_grid`) is shown in the
  [Python model_selection port](../../Python/ml_methods/model_selection/).
- The C++ validation programs (`cpp/ml_methods/validation/`) are
  cross-language checks by design and stay C++-side.

Naming note: the R `standardization/` folder corresponds to the C++
`scaler_methods/` folder.

---

**Last updated**: 2026-07-06
