# Demo 03: T-Rex+GVS on Mixed Blocks

## Purpose

Evaluate T-Rex+GVS on a more realistic grouped design: four contiguous equicorrelated blocks of unequal size, placed in **random order** with **random noise gaps** between and around them, and including one **inactive equicorrelated trap block**. This file merges two benchmark presets that share the same underlying generator, `make_mixed_blocks_dgp()`.

---

## Data Generation Parameters (`make_mixed_blocks_dgp`)

- Four contiguous blocks of sizes $\{20, 50, 80, 65\}$, each following the equicorrelated model
  $$
  X_{ij} = Z_{i,\,g(j)} + \sigma_x\, \xi_{ij}, \qquad \xi_{ij} \sim \mathcal{N}(0,1).
  $$
- **Blocks 0–2 are active** ($\beta_j = 3$); **block 3 is an inactive trap** ($\beta_j = 0$) with the same equicorrelated structure as the active blocks.
- Blocks are placed in a **random order** with **random-width noise gaps** (5 gap slots: before, between, and after the blocks) filled with i.i.d. background columns.
- $n = 200$, $p = 500$, $s = 150$ (sum of the three active blocks).

---

## Two Presets

### Preset 1 — Mixed-Blocks (default $\sigma_x$)
- Part 1: SNR sweep at fixed $\sigma_x = \sqrt{0.01}$ ($\rho \approx 0.99$).
- Part 2: $\rho$ sweep at fixed $\mathrm{SNR}=2.0$, $\rho \in \{0.10,\dots,0.99\}$ (9 levels).
- Part 3: 2-D SNR $\times$ $\rho$ grid.

### Preset 2 — Mixed-Blocks-Rho075 (fixed $\rho = 0.75$)
- Part 1: SNR sweep at fixed $\rho = 0.75$ (derived $\sigma_x$).
- Part 2: $\rho$ sweep at fixed $\mathrm{SNR}=2.0$, with expanded coverage (11 levels, including $0.40$ and $0.60$).
- Part 3: 2-D SNR $\times$ $\rho$ grid. Its $\rho$ axis is $\{0.30, 0.50, 0.70, 0.85, 0.90, 0.99\}$, which replaces Preset 1's $0.95$ with $0.85$ — shifting coverage toward a lower $\rho$ rather than extending it (both presets already reach $0.99$).

Both presets share $\mathrm{SNR} \in \{0.1, 0.2, 0.5, 1.0, 2.0, 5.0\}$ for the SNR sweep and the same 2-D SNR grid $\{0.2, 0.5, 1.0, 2.0, 5.0\}$. Preset 1's 2-D $\rho$ axis is $\{0.30, 0.50, 0.70, 0.90, 0.95, 0.99\}$.

---

## Control Parameters

```
K = 20
tFDR = 0.1
corr_max = 0.98
hc_linkage = Single
lambda2_method = CV_1SE_CCD
MC = 200
```

---

## Methods Compared

- **EN** (`TENET`), **EN+AUG** (`TENET_AUG`), **IEN**.

---

## Output Files

Written to `simulation_results/`:

- `gvs_mixed_blocks_snr.txt` / `.csv`, `gvs_mixed_blocks_rho.txt` / `.csv`, `gvs_mixed_blocks_2d.txt` / `.csv` — Preset 1.
- `gvs_mixed_blocks_rho075_snr.txt` / `.csv`, `gvs_mixed_blocks_rho075_rho.txt` / `.csv`, `gvs_mixed_blocks_rho075_2d.txt` / `.csv` — Preset 2.

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_gvs/demo_trex_gvs_03_mc_sim_mixed_blocks/demo_trex_gvs_03_mc_sim_mixed_blocks
```

Both presets run in the same executable invocation.

---

## Interpretation

- FDR is expected to remain controlled near $\mathrm{tFDR}=0.1$ for **EN**/**EN+AUG** across both presets, since the inactive trap block shares the active blocks' correlation structure and should be excluded correctly once its group is identified as null.
- TPR should increase with SNR as usual, with the random gaps and block order mainly testing that the grouping/clustering logic is invariant to block placement rather than materially changing detection power.
- Comparing Preset 1 (very high default $\rho\approx0.99$) against Preset 2 (moderate $\rho=0.75$) shows how performance degrades as within-block correlation decreases, all else equal.
- As in Demos 01–02, expect **IEN** to be substantially more conservative than the EN-based methods on this design.

---

**Last updated**: 2026-07-08
