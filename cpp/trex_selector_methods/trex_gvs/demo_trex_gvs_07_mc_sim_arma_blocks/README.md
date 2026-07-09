# Demo 07: T-Rex+GVS on Heterogeneous ARMA Blocks

## Purpose

Stress-test T-Rex+GVS when different blocks follow **different time-series structures** rather than a single uniform correlation model. This is the most heterogeneous of the four-block scenarios: each block has its own ARMA specification, testing whether the selector can still recover grouped signal when "the group correlation model" is not one fixed thing.

> **Note**: this demo's `simulation_results/` folder is currently empty in this checkout — the executable has not yet been run to produce committed output. Filenames below are the ones the code writes when executed.

---

## Data Generation Parameters (`make_arma_blocks_dgp`)

Same 4-block geometry as Demo 03/05/06 (sizes $\{20, 50, 80, 65\}$), but each block follows a distinct ARMA process (100-step burn-in per row, columns standardized to unit empirical variance afterward):

- **Block 0** (20 vars, active, $\beta_j=3$): AR(1), $\phi_1 = \texttt{ar\_coef}$.
- **Block 1** (50 vars, active): MA(3), $\theta = (0.5, 0.3, 0.1)$.
- **Block 2** (80 vars, active): ARMA(2,1), $\phi = (\texttt{ar\_coef}, \texttt{ar\_coef}/5)$, $\theta = 0.5$.
- **Block 3** (65 vars, inactive trap): AR(1), $\phi_1 = \texttt{ar\_coef}$.
- Gap columns: i.i.d. $\mathcal{N}(0,1)$.
- $n = 200$, $p = 500$, $s = 150$.

Unlike Demo 06, the sweep parameter here is `ar_coef` (shared by the AR(1)/ARMA(2,1) blocks), **not** $\rho$ — the MA(3) block has no AR component at all.

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

## Three Parts

1. **SNR sweep** at fixed `ar_coef = 0.8`: $\mathrm{SNR} \in \{0.1, 0.2, 0.5, 1.0, 2.0, 5.0\}$.
2. **`ar_coef` sweep** at fixed $\mathrm{SNR}=2.0$: `ar_coef` $\in \{0.10, 0.20, \dots, 0.90, 0.95\}$ (10 levels).
3. **2-D SNR $\times$ `ar_coef` grid**: 5 SNR levels $\times$ 6 `ar_coef` levels.

---

## Output Files (expected)

Written to `simulation_results/` when run:

- `gvs_arma_blocks_snr.txt` / `.csv` — SNR sweep (Part 1).
- `gvs_arma_blocks_arcoef.txt` / `.csv` — `ar_coef` sweep (Part 2).
- `gvs_arma_blocks_2d.txt` / `.csv` — 2-D grid (Part 3).

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_gvs/demo_trex_gvs_07_mc_sim_arma_blocks/demo_trex_gvs_07_mc_sim_arma_blocks
```

---

## Interpretation

- Because the three active blocks have genuinely different correlation structures (AR(1), MA(3), ARMA(2,1)), expect recovery to be **uneven across blocks** — the MA(3) block (short-range, decaying-but-bounded correlation) and the ARMA(2,1) block behave differently from the pure AR(1) block, so aggregate TPR is a blend of these different per-block difficulties.
- FDR should remain reasonably controlled near or below $\mathrm{tFDR}=0.1$, since the trap block (AR(1), same `ar_coef` as the active AR(1) block) is structurally similar to an active block and is the main test of false-group exclusion here.
- Compare against Demo 06 (uniform AR(1) blocks) to see whether heterogeneity across blocks costs additional power beyond what AR(1) decay alone costs.

---

**Last updated**: 2026-07-08
