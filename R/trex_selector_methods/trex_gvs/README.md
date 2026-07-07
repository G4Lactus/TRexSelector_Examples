# T-Rex+GVS: Grouped Variable Selection — R Demonstration Suite

> **Status**: NEW — migrated to the **TRexSelectorNeo** R6 API
> (`TRexGVSSelector$new(..., gvs_control = trex_gvs_control(gvs_type = "EN"|"IEN"))$select()`,
> with `compute_fdp()` / `compute_tpp()`). No longer depends on the CRAN
> `TRexSelector` package.

## Purpose

R demos and Monte Carlo simulations for the **Grouped Variable Selection
T-Rex (T-Rex+GVS)** selector on designs where the active variables cluster
into groups. Each scenario is split into two files:

- `demo_trex_gvs_XX.R` — Part 1: single-run demo with data visualization
  (Plotly correlation heatmaps) and a selector application;
- `simulation_trex_gvs_XX.R` — Parts 2–4: MC simulations comparing **EN vs
  IEN** dummy generation, sweeping SNR (Part 2), a correlation parameter
  (Part 3), and a 2D SNR x correlation phase-transition grid (Part 4).
  Typical fixed settings: n = 200, p = 500, tFDR = 0.1, K = 20, 200 MC runs.

Dependencies: `TRexSelectorNeo`, `plotly`, `parallel`. The scripts source
the shared [../support_generators.R](../support_generators.R) and
[../simulation_utils.R](../simulation_utils.R) (which provides the
`.run_mc_*()` EN/IEN wrappers) plus the local `dgp_*.R` files.

---

## Demos and simulations

| # | Scenario (DGP) | Demo / simulation files | dgp file | cpp counterpart |
|---|---|---|---|---|
| 01 | Hastie (Zou & Hastie 2005): 3 groups of 50 equicorrelated active predictors | [demo_trex_gvs_01.R](demo_trex_gvs_01.R) / [simulation_trex_gvs_01.R](simulation_trex_gvs_01.R) | [dgp_hastie.R](dgp_hastie.R) | `demo_trex_gvs_01_mc_sim_hastie_en_blocks` |
| 02 | Scattered-grouped: 3 groups of 50 randomly scattered across 500 columns | [demo_trex_gvs_02.R](demo_trex_gvs_02.R) / [simulation_trex_gvs_02.R](simulation_trex_gvs_02.R) | [dgp_scattered_grouped.R](dgp_scattered_grouped.R) | `demo_trex_gvs_02_mc_sim_scattered_blocks` |
| 03 | Mixed blocks: 3 active blocks (20/50/80) + 1 inactive trap block (65), random order with noise gaps | [demo_trex_gvs_03.R](demo_trex_gvs_03.R) / [simulation_trex_gvs_03.R](simulation_trex_gvs_03.R) | [dgp_mixed_blocks.R](dgp_mixed_blocks.R) | `demo_trex_gvs_03_mc_sim_mixed_blocks` |
| 04 | Negative traps: active +Z/−Z group plus two inactive +Z/−Z trap groups | [demo_trex_gvs_04.R](demo_trex_gvs_04.R) / [simulation_trex_gvs_04.R](simulation_trex_gvs_04.R) | [dgp_neg_traps.R](dgp_neg_traps.R) | `demo_trex_gvs_04_mc_sim_neg_traps` |
| 05 | Equicorrelated blocks: 4 blocks (20/50/80/65), 3 active + 1 trap, factor-model equicorrelation | [demo_trex_gvs_05.R](demo_trex_gvs_05.R) / [simulation_trex_gvs_05.R](simulation_trex_gvs_05.R) | [dgp_equi_blocks.R](dgp_equi_blocks.R) | — (no direct counterpart; equicorrelated-block scenarios appear in cpp demos 01–04, 08) |
| 06 | Heavy-tailed equi blocks: same 4-block layout with Student-t(3) factors, noise, and shocks | [demo_trex_gvs_06.R](demo_trex_gvs_06.R) / [simulation_trex_gvs_06.R](simulation_trex_gvs_06.R) | [dgp_t3_equi_blocks.R](dgp_t3_equi_blocks.R) | `demo_trex_gvs_05_mc_sim_t3_blocks` |
| 07 | Block-structured AR(1): within-block Cor(X_j, X_k) = rho^\|j−k\| | [demo_trex_gvs_07.R](demo_trex_gvs_07.R) / [simulation_trex_gvs_07.R](simulation_trex_gvs_07.R) | [dgp_ar1_blocks.R](dgp_ar1_blocks.R) | `demo_trex_gvs_06_mc_sim_ar1_blocks` |
| 08 | ARMA mixed-structure blocks: AR(1) / MA(3) / ARMA(2,1) active blocks + AR(1) trap | [demo_trex_gvs_08.R](demo_trex_gvs_08.R) / [simulation_trex_gvs_08.R](simulation_trex_gvs_08.R) | [dgp_arma_blocks.R](dgp_arma_blocks.R) | `demo_trex_gvs_07_mc_sim_arma_blocks` |

The C++ suite additionally contains `demo_trex_gvs_08_mc_sim_block_bench`
(HAC-discovered vs oracle groups on a 100-block benchmark), which is not
mirrored in R.

---

## Additional DGP files

Beyond the eight DGPs used by the current demos, the folder ships further
generators (each provides an interactive `generate_*()` variant and, where
noted, an SNR-calibrated `*_snr()` variant used by the MC wrappers in
[../simulation_utils.R](../simulation_utils.R)):

| File | Description |
|---|---|
| [dgp_gvs.R](dgp_gvs.R) | Collection of GVS benchmark DGPs with a uniform return contract (block equicorrelation and variants) |
| [dgp_unequal_blocks.R](dgp_unequal_blocks.R) | Three contiguous unequal active blocks (20/50/80), rest white noise |
| [dgp_random_blocks.R](dgp_random_blocks.R) | Three active blocks shuffled with random white-noise gaps |
| [dgp_neg_corr.R](dgp_neg_corr.R) | Strict linear negative correlations: active +Z/−Z group and one inactive trap |
| [dgp_factor_model.R](dgp_factor_model.R) | Overlapping latent factor model (K = 4 factors, small active subset) |
| [dgp_hapgen_snr.R](dgp_hapgen_snr.R) | Quasi-hapgen LD-block scenario: 7 irregular LD blocks + weak long-range cross-block LD, p = 500 |

### Known gaps

- `dgp_factor_model.R` has an open TODO: the SNR-calibrated counterpart
  `dgp_factor_model_snr()` is not implemented yet.
- Six DGP files (`dgp_factor_model.R`, `dgp_gvs.R`, `dgp_hapgen_snr.R`,
  `dgp_neg_corr.R`, `dgp_random_blocks.R`, `dgp_unequal_blocks.R`) are not
  sourced by any current demo; they are kept for future scenarios.

---

## Recorded results

Completed simulation runs are documented in
[simulation_results/](simulation_results/):

- [demo_trex_gvs_01_hastie_models.md](simulation_results/demo_trex_gvs_01_hastie_models.md)

---

## Running

```bash
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_01.R        # Part 1
Rscript R/trex_selector_methods/trex_gvs/simulation_trex_gvs_01.R  # Parts 2-4
```

Scripts resolve their own directory, so they run from any working directory.
The simulation scripts expose `run_part_*` flags at the top to enable or
disable individual parts.

---

**Last updated**: 2026-07-08
