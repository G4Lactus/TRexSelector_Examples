# T-Rex+GVS: Grouped Variable Selection Demonstration Suite

## Overview

This folder contains C++ example programs for the **T-Rex Selector with Grouped Variable Selection (T-Rex+GVS)**.

T-Rex+GVS extends the classical T-Rex selector (see [../trex/](../trex/README.md)) to settings where predictors form known or discoverable **groups** of correlated variables — for example genes in the same pathway, or sensor channels driven by a shared latent signal. Instead of voting on individual variables, the selector aggregates dummy-based evidence at the group level, so that a whole group can be selected or rejected together.

The demos in this folder are designed to help users understand:

1. how grouped-signal recovery behaves under different within-group correlation structures (equicorrelated, AR(1), ARMA, heavy-tailed),
2. how group *discovery* via hierarchical agglomerative clustering (HAC) compares to using known ("oracle") groups,
3. how the elastic-net-based (`EN`) and iterative-elastic-net-based (`IEN`) group solvers behave under scattered layouts, correlation traps, and negative correlation structure.

If you are new to this folder, start with **Demo 01** for the canonical equicorrelated-blocks scenario, then continue with **Demo 03** for a more realistic gapped/randomly-ordered block layout, and finish with **Demo 08** for the HAC-vs-oracle group benchmark.

---

## What this folder covers

All demos study sparse variable selection in a Gaussian (or near-Gaussian) linear model where the active variables cluster into groups:

- **equicorrelated blocks** (Demos 01–04, 08): variables within a group share a common latent factor,
- **scattered groups** (Demo 02): group members are permuted across the column space instead of being contiguous,
- **gapped, randomly-ordered blocks with an inactive trap** (Demo 03, and reused by Demos 05–07): four unequal-size blocks, three active and one inactive equicorrelated trap,
- **negative-correlation structure** (Demo 04): groups with sign-flipped halves and spatially separated inactive traps,
- **non-Gaussian / non-equicorrelated structure** (Demos 05–07): Student-t(3) heavy tails, AR(1) Toeplitz covariance, and heterogeneous ARMA processes per block,
- **HAC-discovered vs. oracle groups** (Demo 08): a larger 100-block benchmark comparing automatic group discovery against known group labels, across three within-block truth patterns.

Across demos, the two central questions are: *does the selector recover the correct groups*, and *does the group-level discovery process (HAC clustering) cost much accuracy relative to knowing the true groups in advance*?

---

## Statistical model and notation

We work with the same Gaussian linear model as the classical T-Rex selector,

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(0, \sigma^2 \boldsymbol{I}_n),
$$

but with predictors organized into groups.

### Notation

- $n$: Number of observations, $p$: Number of variables.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$: Design matrix, $\boldsymbol{y} \in \mathbb{R}^n$: Response vector.
- $\boldsymbol{\beta} \in \mathbb{R}^p$: Regression coefficient vector.
- $g(j) \in \{1, \dots, G\}$: Group assignment of variable $j$ (the `prior_groups` vector); background/i.i.d. variables each receive a unique singleton group ID.
- $G$ (or $M$): Number of groups.
- $\mathcal{A} = \{ j : \beta_j \neq 0 \}$: True active-variable support; the corresponding set of **active groups** is $\{ g(j) : j \in \mathcal{A} \}$.
- $\rho$: Within-group correlation (equicorrelated DGPs), induced through a shared latent factor per group.

### Within-group correlation models used across the demos

- **Equicorrelated**: $X_{ij} = \sqrt{\rho}\, Z_{i,g(j)} + \sqrt{1-\rho}\, \varepsilon_{ij}$, with $Z_{\cdot,k} \sim \mathcal{N}(0,1)$ the group-$k$ latent factor.
- **AR(1) Toeplitz** (Demo 06): within-block covariance $\Sigma_k(i,j) = \rho^{|i-j|}$, generated via a Cholesky factor.
- **Heterogeneous ARMA** (Demo 07): each block follows a different ARMA specification (AR(1), MA(3), ARMA(2,1)), swept over an AR coefficient `ar_coef` instead of $\rho$.
- **Heavy-tailed** (Demo 05): the same latent-factor construction as the equicorrelated model, but with Student-t(3)/$\sqrt{3}$ (unit-variance) latent factors and noise.

### Signal-to-noise ratio

As in the classical T-Rex demos, $\mathrm{SNR} = \mathrm{Var}(\boldsymbol{X}\boldsymbol{\beta}) / \mathrm{Var}(\boldsymbol{\varepsilon})$, with $\sigma^2$ calibrated to hit a target SNR for a fixed $(\boldsymbol{X}, \boldsymbol{\beta})$.

---

## T-Rex+GVS specific concepts

Grouped selection introduces a few additional control parameters and types beyond the classical T-Rex selector (from `trex_gvs_simulation_utils.hpp` / `TRexGVSControlParameter`):

- **`GVSType`** — the group-level base method:
  - `EN`: elastic-net-based grouped selection,
  - `IEN`: iterative elastic-net-based grouped selection.
- **`ENSolverType`** (only applies to `GVSType::EN`) — which elastic-net solver drives the dummy experiments:
  - `TENET`: the base Gram-based elastic-net solver,
  - `TENET_AUG`: the augmented-LASSO formulation of the elastic net.
- **`LambdaSelectionMethod`** — how the elastic-net penalty $\lambda_2$ is chosen: `CV_1SE_CCD` (cross-validation, 1-SE rule, via coordinate descent) or `CV_1SE_SVD` (via SVD).
- **`prior_groups`** — the 0-based per-column group-ID vector. When groups are known in advance (an "oracle" run), this vector is taken directly from the data-generating process. When groups must be discovered from the data, HAC clustering (with a `corr_max` correlation threshold and a chosen linkage method) is used to build `prior_groups` automatically — this is what Demo 08 compares directly against oracle groups.

---

## Statistical targets

The same false discovery rate (FDR) and true positive rate (TPR) definitions from the classical T-Rex demos apply here, evaluated either at the **coordinate level** (individual variables) or the **block/group level** (whether an active group was hit, fully recovered, or a null group was falsely activated):

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

Demo 08 additionally reports block-level diagnostics: **block hit rate** (share of active blocks with at least one selected variable), **block FDR** (share of hit blocks that are actually null), **full-block rate** (share of active blocks recovered in their entirety), and a **block-purity** diagnostic (how consistently a true block maps to a single discovered group — always $1.0$ for oracle groups).

---

## Start here

For first use, the following order is recommended:

1. **Demo 01 — Hastie EN Blocks**
   Canonical starting point: three equicorrelated active groups, SNR and $\rho$ sweeps, comparing `EN`, `EN+AUG`, and `IEN`.

2. **Demo 03 — Mixed Blocks**
   A more realistic layout: unequal-size blocks in random order with random gaps, plus one inactive equicorrelated trap block.

3. **Demo 08 — Block Benchmark**
   Compares HAC-discovered groups against oracle groups across three within-block truth patterns (individual / representative / whole-block support).

---

## Demonstration structure

Each demo subfolder typically contains:

- **`demo_*.cpp`**: the source file for the demo.
- **`README.md`**: a description of the scenario and the main interpretation.
- **`simulation_results/`**: generated output files such as text summaries and CSV tables (Demo 08 is console-only and does not write result files).

---

## Demo descriptions

| # | Name | Purpose | Setting | Main parameters |
|---|------|---------|---------|-----------------|
| **01** | Hastie EN Blocks | Canonical equicorrelated grouped-recovery benchmark | 3 active equicorrelated groups (50 vars each) | $n=200$, $p=500$, $s=150$, SNR and $\rho$ sweeps + 2-D grid, EN/EN+AUG/IEN, MC=200 |
| **02** | Scattered Blocks | Tests whether spatial scattering of group members hurts recovery | Same 3 groups as Demo 01, columns randomly permuted | $n=200$, $p=500$, $s=150$, SNR and $\rho$ sweeps + 2-D grid, MC=200 |
| **03** | Mixed Blocks | Realistic gapped, randomly-ordered blocks with an inactive trap | 4 unequal blocks $\{20,50,80,65\}$, 3 active + 1 trap; two presets (default $\sigma_x$ vs. fixed $\rho=0.75$) | $n=200$, $p=500$, $s=150$, SNR and $\rho$ sweeps + 2-D grid, MC=200 |
| **04** | Negative Traps | Recovery of a sign-flipped active group amid spatially separated sign-flipped traps | Active group (0–99) + Trap 1 (100–199) + noise + Trap 2 (300–359) + noise | $n=200$, $p=500$, $s=100$, SNR and $\rho$ sweeps + 2-D grid, MC=200 |
| **05** | T(3) Blocks | Robustness to heavy-tailed, non-Gaussian data | Same 4-block geometry as Demo 03, Student-t(3)/$\sqrt3$ latent factors and noise | $n=200$, $p=500$, $s=150$, SNR and $\rho$ sweeps + 2-D grid, MC=200 |
| **06** | AR(1) Blocks | Robustness to non-equicorrelated (AR(1)) within-block covariance | Same 4-block geometry, $\Sigma_k(i,j)=\rho^{\lvert i-j\rvert}$ | $n=200$, $p=500$, $s=150$, SNR and $\rho$ sweeps + 2-D grid, MC=200 |
| **07** | ARMA Blocks | Robustness to heterogeneous AR/MA/ARMA structure per block | Same 4-block geometry, block-specific AR(1)/MA(3)/ARMA(2,1) processes | $n=200$, $p=500$, $s=150$, SNR and `ar_coef` sweeps + 2-D grid, MC=200 |
| **08** | Block Benchmark | HAC-discovered vs. oracle groups, across truth patterns | 100 equal blocks (size 5), 10 active blocks, 3 truth scenarios (INDIVIDUAL/REPRESENTATIVE/WHOLE) $\times$ 4 methods | $n=200$, $p=500$, $\rho \in \{0.5, 0.9\}$, SNR grid, MC=500 |

---

## Folder contents

```txt
trex_gvs/
  ├── README.md
  ├── CMakeLists.txt
  ├── trex_gvs_dgps.hpp
  ├── trex_gvs_simulation_utils.hpp
  │
  ├── demo_trex_gvs_01_mc_sim_hastie_en_blocks/
  │   ├── demo_trex_gvs_01_mc_sim_hastie_en_blocks.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  ├── demo_trex_gvs_02_mc_sim_scattered_blocks/
  │   ├── demo_trex_gvs_02_mc_sim_scattered_blocks.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  ├── demo_trex_gvs_03_mc_sim_mixed_blocks/
  │   ├── demo_trex_gvs_03_mc_sim_mixed_blocks.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  ├── demo_trex_gvs_04_mc_sim_neg_traps/
  │   ├── demo_trex_gvs_04_mc_sim_neg_traps.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  ├── demo_trex_gvs_05_mc_sim_t3_blocks/
  │   ├── demo_trex_gvs_05_mc_sim_t3_blocks.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  ├── demo_trex_gvs_06_mc_sim_ar1_blocks/
  │   ├── demo_trex_gvs_06_mc_sim_ar1_blocks.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  ├── demo_trex_gvs_07_mc_sim_arma_blocks/
  │   ├── demo_trex_gvs_07_mc_sim_arma_blocks.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  └── demo_trex_gvs_08_mc_sim_block_bench/
      ├── demo_trex_gvs_08_mc_sim_block_bench.cpp
      ├── README.md
      └── simulation_results/   (empty — console-only output)
```

---

## What to expect

Across the equicorrelated, scattered, and trap-block scenarios (Demos 01–04), FDR generally stays near the $\mathrm{tFDR}=0.1$ target while TPR climbs from near zero at low SNR to near-perfect recovery at high SNR, mirroring the classical T-Rex pattern but evaluated at the group level. Robustness scenarios (Demos 05–07) typically show a modest power penalty relative to the plain equicorrelated case as the within-block structure departs from i.i.d. Gaussian equicorrelation (heavy tails, AR(1) decay, heterogeneous ARMA). In the Demo 08 benchmark, HAC-discovered groups (`M1`/`M3`) are expected to track oracle-group performance (`M2`/`M4`) closely when the correlation threshold `corr_max` is well matched to the true within-block correlation, with the largest gaps appearing for the `REPRESENTATIVE` truth scenario where only part of a block is active.

---

## Building the demos

From the C++ workspace root:

```bash
cd TRexSelector_Simulations/cpp
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"
cmake --build build/debug
```

The executables are written to:

```txt
build/debug/bin/
```

---

## Running a demo

For example, to run Demo 01:

```bash
./build/debug/bin/demo_trex_gvs_01_mc_sim_hastie_en_blocks
```

Most demos write both a human-readable `.txt` summary and a tidy `.csv` table into their local `simulation_results/` folder; Demo 08 prints its results to the console only.

---

## Notes for new users

- Start with **Demo 01** to see the baseline equicorrelated-group behavior.
- Use **Demo 03** to see a more realistic, gapped, randomly-ordered block layout with an inactive trap.
- Use **Demos 05–07** if you are interested in robustness to departures from Gaussian equicorrelation (heavy tails, AR(1), ARMA).
- Use **Demo 08** if you are interested in the practical cost of discovering groups via HAC clustering instead of knowing them in advance.
- Check the local `README.md` inside each demo folder for scenario-specific details.

---

## References

- See [../README.md](../README.md) for the category overview of all T-Rex selector variants.
- See [../trex/README.md](../trex/README.md) for the classical (non-grouped) T-Rex selector demos and shared statistical background.
- The R reference implementation and DGP generators for T-Rex+GVS live under `R/trex_selector_methods/trex_gvs/`.
- A diagnostic comparison of HAC stability and EN-solver/scaling interactions is available in `validation/trex_gvs/validation_trex_gvs_01_scaling_solver_comparison.cpp`.

---

**Last updated**: 2026-07-04
