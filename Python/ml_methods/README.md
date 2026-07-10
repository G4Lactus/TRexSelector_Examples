# ML Methods — Python Demos

Data-preprocessing building blocks of the TRexSelector library, exercised
through the `trex_selector_neo.ml_methods` module.

| Subarea | Folder | What it covers | C++ counterpart |
|---|---|---|---|
| HAC clustering | [hac_clustering/](hac_clustering/) | Agglomerative hierarchical clustering with different linkage methods and distance metrics, in-memory and mmap-backed | [cpp/ml_methods/hac_clustering/](../../cpp/ml_methods/hac_clustering/) |
| Standardization | [standardization/](standardization/) | Z-score and Lp-norm column scaling with forward/inverse transforms, in-memory and mmap-backed | [cpp/ml_methods/scaler_methods/](../../cpp/ml_methods/scaler_methods/) |
| PCA | [pca/](pca/) | Principal component analysis: fit, restore round-trip, and out-of-sample transform via the `PCA` class | [cpp/ml_methods/pca/](../../cpp/ml_methods/pca/) |
| SVD | [svd/](svd/) | Truncated singular value decomposition via `SVDSolver` (direct and Gram paths, orthogonality checks) | [cpp/ml_methods/svd/](../../cpp/ml_methods/svd/) |
| Ridge regression | [ridge_regression/](ridge_regression/) | L2-penalized least squares via `RidgeSolver` (primal and dual paths, lambda sweep) | [cpp/ml_methods/ridge_regression/](../../cpp/ml_methods/ridge_regression/) |
| Model selection | [model_selection/](model_selection/) | K-fold cross-validation for ridge (`RidgeCV`) and elastic net (`ElasticNet` path + `ElasticNetCV`), with lambda.min / lambda.1se selection | [cpp/ml_methods/model_selection/](../../cpp/ml_methods/model_selection/) |

---

## Not Yet Ported

Every `ml_methods` demo area from the [cpp/](../../cpp/ml_methods/) tree now
has a Python counterpart. The C++ scaler/HAC validation programs that used to
live under `cpp/ml_methods/validation/` moved to the TRexSelector library test
suite (`TRexSelector/cpp/tests/validation/hac_clustering/`); they are
cross-language reference checks by design, not demos.

---

**Last updated**: 2026-07-06
