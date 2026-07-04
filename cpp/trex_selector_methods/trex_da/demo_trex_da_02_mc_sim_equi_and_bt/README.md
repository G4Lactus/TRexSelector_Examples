# Demo 02: DA-TRex on Equicorrelated and Binary-Tree (BT) Data

## Purpose

Cover two additional dependency structures: **equicorrelation** (`DAMethod::EQUI`) and a **two-level hierarchical block** structure handled via **binary-tree/dendrogram clustering** (`DAMethod::BT`). Corresponds to the R references `demo_trex_da_02_equi.R` and `demo_trex_da_03_bt_equi_blocks.R`.

> **Status**: `simulation_results/` is currently empty — this demo has not yet been run in this checkout.

---

## Data Generation Parameters

### Equicorrelated (`dgp_equi`), Parts 1–2
$$
X_{i,j} = \sqrt{\rho}\, f_i + \sqrt{1-\rho}\, \eta_{i,j}.
$$
$n=300$, $p=1000$, `Random` support.

### Binary-tree / two-level block (`dgp_bt`), Parts 3–4
$$
X_{i,j} = \sqrt{\rho_{\text{between}}}\, f_0 + \sqrt{\rho_{\text{within}}-\rho_{\text{between}}}\, f_{\text{block}(j)} + \sqrt{1-\rho_{\text{within}}}\, \eta_{i,j}.
$$
$n=300$, $p=1000$, $n_{\text{blocks}}=10$, `OnePerBlock` support (one active variable per block).

---

## Four Parts

1. **Equi SNR sweep** — $\rho=0.25$, $\mathrm{SNR}\in\{0.1,0.2,0.5,1.0,2.0,5.0\}$, 3 solvers + base T-Rex comparison.
2. **Equi $\rho$ sweep** — $\mathrm{SNR}=2.0$ fixed, `Random` support.
3. **BT SNR sweep** — $\rho_{\text{within}}=0.5$, $\rho_{\text{between}}=0.1$, `Single` linkage, 3 solvers + base T-Rex comparison.
4. **BT linkage sweep** — $\rho_{\text{within}}=0.8$, $\rho_{\text{between}}=0.1$; outer loop over **Single / Complete / Average** HAC linkage.

All parts: $K=20$, $\mathrm{tFDR}=0.2$, $\mathrm{MC}=200$, `base_seed=2026`.

---

## Output Files (expected)

Written to `simulation_results/` once run — naming pattern `da_trex_mc_da_{equi|bt}_{snr|rho|linkage}.{txt,csv}` (exact stems follow the shared `da_trex_mc_{scenario_tag}` convention; verify exact tag strings in the source if you rerun this demo).

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_da_02_mc_sim_equi_and_bt
```

---

## Interpretation

- Since equicorrelation ties every column to the same shared factor, expect **DA-EQUI** to behave qualitatively like the AR(1) case in Demo 01: reduced FDR vs. base T-Rex, at some TPR cost — but the effect should be more uniform across variables (equicorrelation has no notion of "distance" the way AR(1) does).
- The **BT SNR sweep** (Part 3) is the introduction of the `BT` method used throughout Demos 04–08 — expect its FDR/TPR pattern to depend on how well the HAC dendrogram cut recovers the true block structure.
- The **BT linkage sweep** (Part 4) is the key diagnostic for linkage-method sensitivity: compare Single/Complete/Average directly on the same DGP to see whether one is systematically better suited to the block-equicorrelated structure used here — this same linkage-sensitivity question recurs in Demos 04–07 on more varied block designs.

---

**Last updated**: 2026-07-04
