# Demo 03b: DA+NN Correction Applied to AR(1) Data (Method-Mismatch Study)

## Purpose

A deliberate **misspecification stress test**: apply the **`DAMethod::NN`** (banded/nearest-neighbor) correction to data generated from an **AR(1)** process — a genuinely different (geometrically-decaying, not banded) correlation structure. This studies whether the NN correction, which assumes a finite-range dependency, can still usefully approximate a structure whose correlation never truly vanishes at any distance. Companion to Demo 03 (`dgp_nn`, correctly-specified NN data).

> **Status**: `simulation_results/` is currently empty — this demo has not yet been run in this checkout.

---

## Data Generation Parameters (`dgp_ar1`, same DGP as Demo 01)

$$
X_{i,j} = \rho\, X_{i,j-1} + \sqrt{1-\rho^2}\, \eta_{i,j}.
$$

- $n=300$, $p=1000$, $\rho=0.7$ — identical DGP parameters to Demo 01, but here the selector is configured with `da_ctrl.method = DAMethod::NN` instead of `AR1`.

---

## Three Parts

1. **SNR sweep** — $\rho=0.7$ fixed; `CappedSpread(max_gap=20)` and `Random` support; 3 solvers + base T-Rex comparison.
2. **$\rho$ sweep** — $\mathrm{SNR}=2.0$ fixed; dual support.
3. **2D SNR × $\rho$ sweep** — all 3 solvers; `CappedSpread` and `Random` sub-sections.

All parts: $s=10$, $K=20$, $\mathrm{tFDR}=0.2$, $\mathrm{MC}=200$, `base_seed=2026`.

---

## Output Files (expected)

Written to `simulation_results/` once run, following the shared `da_trex_mc_{scenario_tag}.{txt,csv}` naming convention (verify exact tag strings in the source if you rerun this demo).

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_da_03b_mc_sim_nn_ar
```

---

## Interpretation

- This demo directly tests **robustness to model misspecification**: NN assumes correlation vanishes beyond $\kappa$ columns, while AR(1) correlation decays but never truly reaches zero. Expect DA-NN to recover *some* of the FDR-control benefit seen with correctly-specified `AR1` correction in Demo 01, but likely **less completely** — some residual correlation beyond the NN correction's effective window should leak through as extra false discoveries relative to Demo 01's matched `AR1`-on-`AR1` results.
- Read this demo side-by-side with **Demo 01** (same DGP, correctly-specified `AR1` correction) and **Demo 03** (correctly-specified `NN` correction, different DGP) to triangulate how much of DA-TRex's benefit comes from having the *exactly right* dependency model vs. an *approximately reasonable* one.
- If DA-NN performs nearly as well as DA-AR1 here, that would suggest the method is reasonably robust to this kind of correlation-model mismatch; a large gap would suggest the correction is fairly sensitive to matching the true dependency structure.

---

**Last updated**: 2026-07-04
