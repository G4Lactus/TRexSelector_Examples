# DA-TRex: Dependency-Aware T-Rex Demonstration Suite

## Overview

This folder contains C++ example programs for **DA-TRex**, the T-Rex Selector extended with **dependency-aware (DA) correction** for correlated/dependent design matrices. Where the classical T-Rex selector (see [../trex/](../trex/README.md)) assumes exchangeable, weakly correlated predictors, DA-TRex explicitly models the correlation structure among nearby or grouped variables so that a single true signal does not get "smeared" into false discoveries among its correlated neighbors.

The demos cover four dependency structures and their matching DA methods: **AR(1)** (autoregressive/Toeplitz decay), **EQUI** (equicorrelation), **NN** (banded/nearest-neighbor), and **BT** (binary-tree/dendrogram-based hierarchical block detection via HAC clustering).

---

## What this folder covers

- **AR(1) Toeplitz correlation** (Demo 01): the classical case DA-TRex was designed for — correlation decays geometrically with column distance.
- **Equicorrelation and hierarchical blocks** (Demo 02): compound-symmetry correlation, and a two-level block structure handled via BT.
- **Banded/nearest-neighbor correlation** (Demo 03) and a **method-mismatch stress test** (Demo 03b: applying NN correction to AR(1) data).
- **Block-diagonal AR(1) designs**, with and without appended white-noise columns (Demos 04–05, 07), using BT aggregation with a sweep over HAC linkage methods.
- **Heavy-tailed (Student-t) block designs** (Demos 06–07): testing DA-TRex robustness when predictors and/or noise depart from Gaussianity.
- **Multi-level nested group structures** (Demo 08): three-level hierarchical latent factors with a non-exchangeable Toeplitz leaf layer.

A correctness diagnostic comparing the BT clustering pipeline step-by-step against an R reference implementation lives in [../validation/trex_da/](../validation/trex_da/README.md), not in this demo suite.

---

## Statistical model and notation

DA-TRex uses the same Gaussian linear model as the classical T-Rex selector,

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(0, \sigma^2 \boldsymbol{I}_n),
$$

but $\boldsymbol{X}$'s columns are no longer (close to) exchangeable — they follow one of several explicit dependency structures:

- **AR(1) / Toeplitz**: $\Sigma_{jk} = \rho^{|j-k|}$, built via the recursion $X_{i,j} = \rho X_{i,j-1} + \sqrt{1-\rho^2}\,\eta_{i,j}$.
- **Equicorrelated**: $X_{i,j} = \sqrt{\rho}\, f_i + \sqrt{1-\rho}\,\eta_{i,j}$, a single shared latent factor $f_i$ across all columns.
- **Banded / MA($\kappa$))**: $X_{i,j} = \sum_{l=0}^{\kappa} \theta_l\, \eta_{i,j+l}$, with $\theta_l \propto \rho^l$ normalized to unit variance — only the $\kappa$ nearest neighboring columns are correlated.
- **Binary-tree / two-level block**: $X_{i,j} = \sqrt{\rho_{\text{between}}}\, f_0 + \sqrt{\rho_{\text{within}} - \rho_{\text{between}}}\, f_{\text{block}(j)} + \sqrt{1-\rho_{\text{within}}}\,\eta_{i,j}$, requiring $\rho_{\text{within}} > \rho_{\text{between}}$.
- **Multi-level nested groups**: $X_{i,j} = \sum_{x} \sqrt{\rho_x - \rho_{x+1}}\, f_{g_x(j)} + \sqrt{1-\rho_0}\,\eta_{i,j}$, with strictly decreasing $\rho$ levels from fine to coarse groupings.

Signal-to-noise ratio follows the same convention as elsewhere: $\mathrm{SNR} = \mathrm{Var}(\boldsymbol{X}\boldsymbol{\beta})/\mathrm{Var}(\boldsymbol{\varepsilon})$.

---

## DA-TRex specific concepts

- **`DAMethod`** — the dependency-correction strategy: `AR1`, `EQUI`, `NN`, `BT`. Each targets the matching correlation structure above; DA-TRex uses this structure to widen or narrow the "correction window" around each candidate variable when accumulating dummy-based voting evidence.
- **BT linkage method** (`hac::LinkageMethod`) — for the `BT` method, variables are grouped via hierarchical agglomerative clustering; several demos sweep **Single**, **Complete**, and **Average** linkage to study sensitivity to this choice.
- **Support placement policies** (`SupportPolicy`): `CappedSpread` (deterministic, evenly spaced with a capped gap — used to probe how spacing interacts with a correction window), `Random` (redrawn per trial), `Clustered` (stochastic small groups), `OnePerBlock` (one active variable per block — the natural choice for block-structured DGPs).
- **Base T-Rex comparison row**: several demos additionally run the *classical* T-Rex selector (no DA correction) side by side with DA-TRex on the same correlated data, to directly show the FDR-control benefit (or TPR cost) of the dependency-aware correction.

---

## Statistical targets

Same FDR/TPR definitions as throughout this project:

$$
\mathrm{FDR}
=
\mathbb{E}
\left[
\frac{|\widehat{\mathcal{A}} \setminus \mathcal{A}|}
{\max\{1, |\widehat{\mathcal{A}}|\}}
\right],
\qquad
\mathrm{TPR}
=
\mathbb{E}
\left[
\frac{|\mathcal{A} \cap \widehat{\mathcal{A}}|}
{\max\{1, |\mathcal{A}|\}}
\right].
$$

See [../validation/trex_da/](../validation/trex_da/README.md) for a **correctness diagnostic** — not a statistical-power study — that compares C++ and R dendrogram heights, group memberships, and selection results on one identical dataset to check they agree, rather than measuring FDR/TPR over many trials.

---

## Start here

1. **Demo 01 — AR(1)**: the foundational scenario, including a direct DA-TRex vs. base-T-Rex comparison and a "gap × rho" sweep showing when correlated neighbors start hurting the classical selector.
2. **Demo 04 — BT on block AR(1)**: introduces the `BT` method and HAC linkage sweep on a clean block-diagonal design.
3. **Demo 08 — Prior groups**: the most structurally rich DGP (three nested levels of correlation).
4. **[../validation/trex_da/](../validation/trex_da/README.md)** — if you want to verify the BT clustering pipeline itself against an R reference on one fixed dataset.

---

## Demo descriptions

| # | Name | DA Method | DGP | Key parameters | Status |
|---|------|-----------|-----|-----------------|--------|
| **01** | AR(1) | `AR1` | `dgp_ar1` | $n=300,p=1000,s=10,\rho=0.7$; SNR & $\rho$ sweeps + 2D gap×$\rho$ | Real data |
| **02** | Equi + BT | `EQUI`, `BT` | `dgp_equi`, `dgp_bt` | $n=300,p=1000$; equi $\rho=0.25$; BT $n_{\text{blocks}}=10,\rho_w=0.5,\rho_b=0.1$ | Not yet run |
| **03** | NN | `NN` | `dgp_nn` | $n=300,p=1000,\kappa=3,\rho=0.7$; incl. 2D $\kappa\times\rho$ | Not yet run |
| **03b** | NN (mismatch test) | `NN` | `dgp_ar1` | Same as 01 but corrected with NN instead of AR1 — misspecification stress test | Not yet run |
| **04** | BT | `BT` | `dgp_ar1_block` | $n=150,M=5,Q=5(p=25,s=5),\rho=0.7$; linkage×{SNR,$\rho$,Q,M} sweeps | Real data |
| **05** | BT + white noise | `BT` | `dgp_ar1_block_white` | $n=300,p_{\text{total}}=1000,M=5,Q=5(p_{ar}=25),\rho=0.7$; linkage×{SNR,$\rho$,Q,M} sweeps | Real data |
| **06** | BT, heavy-tailed | `BT` | `dgp_block_toeplitz_hvt` | $n=150,M=5,Q=5,\rho=0.8,\nu=3$; Gauss/Heavy × linkage × {SNR,$\rho$,Q,M,tFDR} | Real data |
| **07** | BT, heavy-tailed + white | `BT` | `dgp_ht_block_white` | $n=150,p_{\text{total}}=500,M=5,Q=5,\rho=0.8,\nu=3$; Gauss/Heavy × linkage × {SNR,$\rho$,Q,M,tFDR} | Real data |
| **08** | Prior groups | `BT`-style prior groups | `dgp_groups_toeplitz_leaf` | $n=300,p=1000,s=10$; 3-level groups $\{10,50,250\}$, $\rho=\{0.55,0.25,0.10\}$, $\phi=0.5$ | Real data |

(A correctness diagnostic previously numbered Demo 09 now lives in [../validation/trex_da/validation_trex_da_01_bt_dendro_diag/](../validation/trex_da/validation_trex_da_01_bt_dendro_diag/README.md).)

---

## Folder contents

```txt
trex_da/
  ├── README.md
  ├── CMakeLists.txt
  ├── dgp_generators.hpp
  ├── trex_da_sim_common.hpp
  ├── demo_trex_da_01_mc_sim_ar1/
  ├── demo_trex_da_02_mc_sim_equi_and_bt/
  ├── demo_trex_da_03_mc_sim_nn/
  ├── demo_trex_da_03b_mc_sim_nn_ar/
  ├── demo_trex_da_04_mc_sim_bt_ar1_block/
  ├── demo_trex_da_05_mc_sim_bt_ar1_block_sweeps/
  ├── demo_trex_da_06_mc_sim_bt_ht_block_sweeps/
  ├── demo_trex_da_07_mc_sim_bt_ht_block_white/
  └── demo_trex_da_08_mc_sim_groups/
```

(Each demo subfolder above also contains its own `.cpp` file, `README.md`, and `simulation_results/`.)

---

## What to expect

For demos with real committed data (01, 04–08), a consistent pattern emerges: DA-TRex's FDR stays much closer to well-controlled than the classical (no-DA) T-Rex selector on the same correlated data — for example, in Demo 01 at $\mathrm{SNR}=2.0$, DA-TRex achieves $\mathrm{FDR}\approx0.059$ vs. base T-Rex's $\mathrm{FDR}\approx0.29$ at the same TPR ballpark — at some cost in TPR, especially at low SNR. Heavy-tailed and grouped-DGP demos (06–08) show correspondingly higher realized FDR relative to the $\mathrm{tFDR}$ target, reflecting the added difficulty of heavy tails and multi-level dependency. For demos without committed output yet (02, 03, 03b), treat any statements as expectations to verify once run, not confirmed results.

---

## Building the demos

```bash
cd TRexSelector_Examples/cpp
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"
cmake --build build/debug
```

---

## Running a demo

```bash
./build/debug/bin/demo_trex_da_01_mc_sim_ar1
```

Most demos write both `.txt` and `.csv` files per scenario into their local `simulation_results/` folder, following the naming pattern `da_trex_mc_{scenario_tag}.{txt,csv}`.

---

## Notes for new users

- Start with **Demo 01** to see the core AR(1) scenario and the DA-vs-base-T-Rex comparison.
- Use **Demos 04–07** to explore the `BT` method and HAC linkage sensitivity on block-structured designs, including heavy-tailed robustness.
- Use **Demo 08** for the richest hierarchical dependency structure.
- See [../validation/trex_da/](../validation/trex_da/README.md) if you need to debug or verify the BT clustering pipeline itself against the R reference.
- Cross-check against the R reference implementation in `R/trex_selector_methods/trex_da/` (extensive: per-demo `.R` scripts plus shared `dgp_*.R`, `simulation_utils.R`, `support_generators.R`).

---

## References

- See [../README.md](../README.md) for the category overview of all T-Rex selector variants.
- See [../trex/README.md](../trex/README.md) for the classical (non-DA) T-Rex selector and shared statistical background.

---

**Last updated**: 2026-07-04
