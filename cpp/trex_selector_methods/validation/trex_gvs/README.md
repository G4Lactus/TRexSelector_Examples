# Validation: T-Rex+GVS Scaling / Solver Comparison

## Overview

This folder holds correctness diagnostics for the **T-Rex+GVS** (grouped variable selection) implementation. It currently contains a single validation program.

---

## `validation_trex_gvs_01_scaling_solver_comparison`

Diagnoses a GVS / elastic-net (EN) scaling interaction by separating a possible
hierarchical-clustering failure from a solver/scaling effect in the EN branch.
It runs two parts on the Hastie equicorrelated-blocks DGP (`make_hastie_dgp`,
$n=200$, $p=500$, $\rho\approx0.99$):

- **Part 1 — clustering integrity**: checks whether the HAC group labels are
  unchanged under `ScalingMode::L2` versus `ScalingMode::ZSCORE` on one fixed
  dataset.
- **Part 2 — scaling × solver grid**: a Monte Carlo 2×2 grid over
  `{L2, ZSCORE} × {TENET, TENET_AUG}`, reporting mean FDR, mean TPR,
  `lambda2_used`, `T_stop`, `num_dummies`, and `M_found`.

The console banner also prints the R reference values for the matched setting
(EN mean_FDP $\approx 0.1174$, TPP $\approx 1.000$) for comparison.

---

## Running

```bash
./build/debug/bin/trex_selector_methods/validation/trex_gvs/validation_trex_gvs_01_scaling_solver_comparison
```

This program prints its results to the console; it does not require external
reference-data files.

---

## References

- The T-Rex+GVS demo suite lives in [../../trex_gvs/](../../trex_gvs/README.md).
- DGP generators and the MC harness are shared from
  `../../trex_gvs/trex_gvs_dgps.hpp` and
  `../../trex_gvs/trex_gvs_simulation_utils.hpp`.

---

**Last updated**: 2026-07-08
