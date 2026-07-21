# Scaler Methods: Demonstration Suite

## Overview

This folder contains Python examples for **column-wise scaling and normalization** utilities in the
`trex_selector_neo.ml_methods` module.

These preprocessing methods are used to standardize features before downstream tasks such as penalized regression, and
other numerical linear-algebra workflows. The current demo illustrates both variance-based standardization and
norm-based normalization. The C++ counterparts live in
[cpp/ml_methods/scaler_methods/](../../../cpp/ml_methods/scaler_methods/).

---

## Key concepts

Given a data matrix $\boldsymbol{X} \in \mathbb{R}^{n \times p}$, scaling is applied column by column.

Two core transformations are demonstrated:

### Z-score scaling

$$
\tilde{x}_{ij} = \frac{x_{ij} - \bar{x}_j}{s_j},
$$

where $\bar{x}_j$ is the column mean and $s_j$ is the column scale.

### Lp-norm scaling

$$
\tilde{\boldsymbol{x}}_j =
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_p},
$$

when centering is enabled.

---

## Start here

| Folder | Purpose |
| ------ | ------- |
| [demo_mlm_scaler_01/](demo_mlm_scaler_01/README.md) | Demonstrates `ZScoreScaler` and `LpNormScaler`, including inverse transforms, serialization, and memory-mapped usage |

---

## Running

```bash
python demo_mlm_scaler_01/demo_mlm_scaler_01.py
```

---

## History note

This folder was previously named `standardization/` and contained a separate mmap-only script
(`demo_standardization_02_mmap.py`). The C++ demo covers memory-mapped scaling as Subproblem 4 of
`demo_mlm_scaler_01`, so the mmap walkthrough was folded into that demo and the standalone mmap variant was retired.

---

**Last updated**: 2026-07-21
