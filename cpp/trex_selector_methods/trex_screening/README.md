# Screen-TRex: Variable Screening Demonstration Suite

## Overview

This folder contains C++ example programs for **Screen-TRex**, the T-Rex Selector extended with a variable-**screening** step for ultra-high-dimensional data. Screen-TRex is designed for settings where $p$ is so large that running the classical T-Rex selector directly is impractical (e.g. genome-wide or biobank-scale designs), by first thresholding a dummy-based voting statistic to screen down to a smaller candidate set.

The demos in this folder are designed to help users understand:

1. how the basic **Ordinary** (majority-vote) and **Bootstrap-CI** screening thresholds work, in-memory and memory-mapped,
2. how **dependency-aware (DA)** screening variants recover power under correlated designs (AR(1), equicorrelated, block-equicorrelated) where naive screening degrades,
3. how the **Biobank** workflow (`BiobankScreenTRex`, "Algorithm 1") adaptively routes each of many phenotypes through Ordinary → Bootstrap-CI → classical T-Rex fallback,
4. how the choice of underlying T-Rex solver (TLARS, TAFS, TOMP) interacts with the screening method.

> **Status note**: every demo in this folder currently has an **empty `simulation_results/` folder** — none of the six executables have been run and committed yet in this checkout. All statements below describe the intended/expected behavior of each demo based on its code and control parameters, not confirmed empirical results.

If you are new to this folder, start with **Demo 01** for the basic Ordinary/Bootstrap comparison, then **Demo 03** to see why dependency-aware screening matters under correlation, then **Demo 04** for the multi-phenotype Biobank workflow.

---

## What this folder covers

- **Basic screening** (Demos 01–02): Ordinary vs. Bootstrap-CI thresholding, in-memory and memory-mapped.
- **Correlation robustness** (Demo 03): naive ("Ordinary") screening vs. dependency-aware (DA) variants under AR(1), equicorrelated, and block-equicorrelated designs.
- **Multi-phenotype biobank screening** (Demos 04–05): Algorithm 1's adaptive Ordinary → Bootstrap-CI → T-Rex-fallback routing, applied to several phenotypes sharing one design matrix, in-memory and memory-mapped.
- **Solver sensitivity** (Demo 06): whether the underlying T-Rex solver backend (TLARS, TAFS, TOMP) changes screening performance.

---

## Statistical model and notation

Screen-TRex operates on the same Gaussian linear model as the classical T-Rex selector,

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(0, \sigma^2 \boldsymbol{I}_n),
$$

but is designed for the regime where $p$ is very large relative to what a direct T-Rex run can handle, so a screening step is applied first to reduce the candidate variable set.

### Notation

- $n$, $p$: observations and variables; $\mathcal{A} = \{j : \beta_j \neq 0\}$: true active support.
- $\Phi_j$: the dummy-based voting proportion for variable $j$ (the same evidence T-Rex accumulates internally); the **Ordinary** screening rule selects $\{ j : \Phi_j > 0.5 \}$.
- $\rho$: pairwise correlation strength in the correlated-design demos (Demo 03), used either as an AR(1) lag-1 correlation, an equicorrelation, or a within-block equicorrelation.

### Signal-to-noise ratio

As elsewhere, $\mathrm{SNR} = \mathrm{Var}(\boldsymbol{X}\boldsymbol{\beta}) / \mathrm{Var}(\boldsymbol{\varepsilon})$, calibrated via $\sigma$ for a fixed $(\boldsymbol{X}, \boldsymbol{\beta})$.

---

## Screen-TRex specific concepts

- **`ScreenTRexMethod`** — the screening variant: `TREX` (ordinary voting-based screening, used with or without bootstrap), and the dependency-aware variants `TREX_DA_AR1`, `TREX_DA_EQUI`, `TREX_DA_BLOCK_EQUI` (Demo 03 only).
- **Ordinary vs. Bootstrap-CI** (`ScreenTRexControlParameter::use_bootstrap_CI`): the **Ordinary** rule thresholds the voting proportion at $\Phi_j > 0.5$; the **Bootstrap-CI** rule instead builds a bootstrap confidence band (`R_boot = 1000` replicates, grid granularity `ci_grid_step = 0.001`) around the estimated FDR curve and picks a threshold from that band.
- **Dependency-aware (DA) screening** (`rho_thr_DA = 0.02`, `n_blocks`): pre-groups strongly correlated variables (above the `rho_thr_DA` correlation threshold, or into `n_blocks` blocks for the block variant) before voting, so that correlated redundant variables do not split evidence and suppress power the way they would under Ordinary screening.
- **Biobank / "Algorithm 1"** (`BiobankScreenTRexControl`, Demos 04–05): screens **multiple phenotypes against one shared design matrix $\boldsymbol{X}$**. For each phenotype, it tries **Screen-TRex Ordinary** first; if the estimated FDR falls outside the acceptable window (`lower_bound_FDR = 0.05`, `upper_bound_FDR = 0.15`), it tries **Screen-TRex Bootstrap-CI**; if that still fails, it falls back to the **classical T-Rex selector** at a fixed `target_FDR_trex = 0.10`. The fraction of phenotypes/trials routed to each method is reported as **"Usage %"**.

---

## Statistical targets

The same FDR/TPR definitions used throughout this project apply here:

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

Screen-TRex additionally reports an **"Estimated FDR"** (the screening procedure's own internal FDR estimate) alongside the realized FDR/TPR, so users can check whether the screening threshold's self-reported FDR tracks the actual empirical FDR. In the biobank demos, FDP/TPP are accumulated **conditionally** on which method a phenotype was routed to, alongside the unconditional **"Usage %"**.

---

## Start here

1. **Demo 01 — Screen-TRex (in-memory)**: basic Ordinary vs. Bootstrap-CI comparison.
2. **Demo 03 — Screen-TRex on correlated designs**: see why dependency-aware screening matters.
3. **Demo 04 — Biobank Screen-TRex (in-memory)**: multi-phenotype adaptive routing.

---

## Demo descriptions

| # | Name | Purpose | Setting | Main parameters |
|---|------|---------|---------|-----------------|
| **01** | Screen-TRex | Baseline Ordinary vs. Bootstrap-CI screening | In-memory, $n=300$, $p=1000$ | $s=5$ (single run) / $s=10$ (MC), SNR sweep 8 points, MC=500 |
| **02** | Screen-TRex (mmap) | Same as Demo 01 with memory-mapped dummy storage | Memory-mapped, $n=300$, $p=1000$ | Same as Demo 01 |
| **03** | Screen-TRex Correlated | Ordinary vs. DA variants under AR(1)/equicorrelated/block-equicorrelated designs | $n=300$, $p=1000$, $s=10$ | 9 parts: 4 single-runs + SNR sweeps + $\rho$ sweeps, MC=500 |
| **04** | Biobank Screen-TRex (in-memory) | Multi-phenotype adaptive routing (Algorithm 1) | $n=300$, $p=1000$, shared $\boldsymbol{X}$ | 1 or 5 phenotypes; MC SNR sweeps, MC=500 |
| **05** | Biobank Screen-TRex (mmap) | Same as Demo 04 with memory-mapped dummy storage | Memory-mapped | Same as Demo 04 |
| **06** | Screen-TRex Solvers | Solver backend comparison (TLARS/TAFS/TOMP) | $n=300$, $p=1000$, $s=10$ | 8-point SNR sweep, MC=500, 3 or 6 method series |

---

## Folder contents

```txt
trex_screening/
  ├── README.md
  ├── CMakeLists.txt
  ├── demo_trex_scr_common.hpp
  ├── demo_trex_scr_bb_common.hpp
  │
  ├── demo_trex_scr_01_screen_trex/
  │   ├── demo_trex_scr_01_screen_trex.cpp
  │   ├── README.md
  │   └── simulation_results/   (empty — not yet run)
  │
  ├── demo_trex_scr_02_screen_trex_mmap/
  │   ├── demo_trex_scr_02_screen_trex_mmap.cpp
  │   ├── README.md
  │   └── simulation_results/   (empty — not yet run)
  │
  ├── demo_trex_scr_03_screen_trex_correlated/
  │   ├── demo_trex_scr_03_screen_trex_correlated.cpp
  │   ├── README.md
  │   └── simulation_results/   (empty — not yet run)
  │
  ├── demo_trex_scr_04_biobank_screen_trex_inmem/
  │   ├── demo_trex_scr_04_biobank_screen_trex_inmem.cpp
  │   ├── README.md
  │   └── simulation_results/   (empty — not yet run)
  │
  ├── demo_trex_scr_05_biobank_screen_trex_mmap/
  │   ├── demo_trex_scr_05_biobank_screen_trex_mmap.cpp
  │   ├── README.md
  │   └── simulation_results/   (empty — not yet run)
  │
  └── demo_trex_scr_06_screen_trex_solvers/
      ├── demo_trex_scr_06_screen_trex_solvers.cpp
      ├── README.md
      └── simulation_results/   (empty — not yet run)
```

---

## What to expect

Since no demo has been executed yet in this checkout, treat the following as expectations to verify, not confirmed results: Ordinary screening should control FDR reasonably well on uncorrelated or weakly correlated designs but is expected to lose power (and potentially inflate estimated vs. realized FDR) as within-predictor correlation grows, which is exactly what the DA variants in Demo 03 are designed to correct. In the biobank demos, Ordinary screening should dominate "Usage %" at high SNR, with Bootstrap-CI and the T-Rex fallback invoked more often as SNR drops or per-phenotype signal becomes harder to separate from noise. Demo 06 should reveal whether solver choice mainly affects speed (as in the classical T-Rex demos) or also shifts the FDR/TPR tradeoff under screening.

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
./build/debug/bin/trex_selector_methods/trex_screening/demo_trex_scr_01_screen_trex/demo_trex_scr_01_screen_trex
```

(The build mirrors the source-tree layout under `bin/`, so each demo executable lives at `bin/trex_selector_methods/trex_screening/<demo_folder>/<demo_name>`.)

Each demo writes both a `.txt` summary and a tidy `.csv` table into its local `simulation_results/` folder once run.

---

## Notes for new users

- Start with **Demo 01** to see the baseline Ordinary/Bootstrap-CI behavior.
- Use **Demo 03** if you are interested in correlation robustness and dependency-aware screening.
- Use **Demos 04–05** for the multi-phenotype biobank workflow and its adaptive method routing.
- Use **Demo 06** if you are interested in solver-choice sensitivity.
- No R reference implementation or validation suite exists yet for Screen-TRex (unlike `trex/` and `trex_gvs/`); this is worth keeping in mind when cross-checking results against other languages.

---

## References

- See [../README.md](../README.md) for the category overview of all T-Rex selector variants.
- See [../trex/README.md](../trex/README.md) for the classical T-Rex selector and shared statistical background.

---

**Last updated**: 2026-07-08
