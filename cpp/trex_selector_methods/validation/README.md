# Validation: T-Rex Selector Methods Correctness Checks

## Overview

This folder contains **validation programs** for the `trex_selector_methods` module.

Unlike lightweight unit tests, these validations are designed as heavier numerical equivalence and invariance checks. They typically rerun full selector workflows on fixed or reference datasets and compare the C++ outputs against trusted reference implementations, especially R-based baselines (`TRexSelector`).

The goals of this folder are:

1. to verify numerical correctness of the C++ selector implementations,
2. to confirm cross-language agreement with reference R code,
3. to test full workflow stages (e.g. clustering, scaling, solver choice) that are too heavy or too data-dependent for ordinary unit tests.

---

## Current categories

| Category | Folder | Purpose |
|----------|--------|---------|
| Classical T-Rex | `trex/` *(README pending)* | Invariance checks for the classical T-Rex selector (e.g. column scaling) |
| Dependency-Aware T-Rex (DA-TRex) | [trex_da/](trex_da/README.md) | Cross-checks the DA-BT hierarchical clustering pipeline against an R reference implementation |
| Grouped Variable Selection (T-Rex+GVS) | [trex_gvs/](trex_gvs/README.md) | Invariance checks for GVS-based grouping and solver comparisons |
| Sparse PCA (T-Rex SPCA) | `trex_spca/` *(README pending)* | Solver comparisons, scaling invariance, and R-dump pipeline cross-checks for T-Rex SPCA |

---

## Output convention

Validation programs write any optional generated files to `validation_results/` rather than `simulation_results/`.

---

## Running

Where a subfolder provides a README (see [trex_da/](trex_da/README.md) and [trex_gvs/](trex_gvs/README.md)), it lists the exact executable and any required reference-data location. Subfolders still marked *(README pending)* above do not yet have one.

---

**Last updated**: 2026-07-08
