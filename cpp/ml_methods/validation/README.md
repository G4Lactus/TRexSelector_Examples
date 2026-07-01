# Validation: ML Methods Correctness Checks

## Overview

This folder contains **validation programs** for the `ml_methods` module.

Unlike lightweight unit tests, these validations are designed as heavier numerical equivalence checks. They typically rerun full algorithmic workflows on fixed datasets and compare the C++ outputs against trusted reference implementations, especially R-based baselines.

The goals of this folder are:

1. to verify numerical correctness of the C++ implementations,
2. to confirm cross-language agreement with reference R code,
3. to test full workflow stages that are too heavy or too data-dependent for ordinary unit tests.

---

## Current categories

| Category | Folder | Purpose |
|----------|--------|---------|
| HAC Clustering | [hac_clustering/](hac_clustering/README.md) | Cross-checks hierarchical clustering and GVS dummy-generation logic against R reference outputs |

---

## Output convention

Validation programs write any optional generated files to `validation_results/` rather than `simulation_results/`.

---

## Running

See the README inside each validation subfolder for the exact executable and required reference-data location.

---

**Last updated**: 2026-07-01
