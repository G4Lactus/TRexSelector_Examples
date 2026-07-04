# ML Methods: Demos and Validation

## Overview

This folder contains the `trex::ml_methods` layer of the C++ project.

It provides the numerical and preprocessing components that support higher-level T-Rex workflows, including column scaling, singular value decomposition, principal component analysis, ridge regression, model selection, hierarchical agglomerative clustering, and validation against R reference implementations.

Each method family has its own subfolder. Demo programs illustrate usage and interface behavior, while validation programs perform heavier correctness checks on fixed reference datasets.

---

## Category overview

| Category | Folder | Purpose |
|----------|--------|---------|
| **Scaler Methods** | [scaler_methods/](scaler_methods/README.md) | Column normalization utilities, including Z-score and Lp-norm scaling |
| **PCA** | [pca/](pca/README.md) | Principal component analysis |
| **SVD** | [svd/](svd/README.md) | Singular value decomposition for stable matrix computations |
| **Ridge Regression** | [ridge_regression/](ridge_regression/README.md) | L2-penalized least-squares regression |
| **Model Selection** | [model_selection/](model_selection/README.md) | Cross-validation workflows for regularized regression methods |
| **HAC Clustering** | [hac_clustering/](hac_clustering/README.md) | Hierarchical agglomerative clustering, including large-scale memory-mapped workflows |
| **Validation** | [validation/](validation/README.md) | Heavy numerical correctness checks against R reference outputs |

---

## Architecture context

These components form the `ml_methods` layer in the project structure:

```txt
utils -> ml_methods -> tsolvers -> trex_selector_methods -> demos
```

In particular:

- scaler methods support feature normalization before downstream modeling,
- SVD and ridge-related tools support stable regression and dimensionality-reduction workflows,
- HAC clustering supports grouped-variable workflows such as GVS,
- validation programs check C++ behavior against trusted R reference pipelines.

---

## Building

```bash
cd TRexSelector_Examples/cpp
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"
cmake --build build/debug
```

---

## Running

Run any demo or validation executable from the build tree. For example:

```bash
./build/debug/bin/ml_methods/pca/demo_mlm_pca_01/demo_mlm_pca_01
```

For validation executables, see the corresponding subfolder README for required reference-data inputs.

---

## Folder contents

```txt
ml_methods/
  ├── README.md
  ├── CMakeLists.txt
  ├── scaler_methods/
  ├── pca/
  ├── svd/
  ├── ridge_regression/
  ├── model_selection/
  ├── hac_clustering/
  └── validation/
```

---

**Last updated**: 2026-07-01
