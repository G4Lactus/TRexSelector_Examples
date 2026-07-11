# DA-TRex: Dependency-Aware T-Rex — R Demonstration Suite

> **Status**: NEW — migrated to the **TRexSelectorNeo** R6 API
> (`TRexDASelector$new(..., da_control = trex_da_control(...))$select()`, with
> `compute_fdp()` / `compute_tpp()`). No longer depends on the CRAN
> `TRexSelector` package.

## Purpose

R demos and Monte Carlo simulations for the **Dependency-Aware T-Rex
(DA-TRex)** selector across correlated designs: AR(1) Toeplitz,
equicorrelation, nearest-neighbor (banded MA), block-diagonal AR(1) (with and
without white-noise augmentation), and heavy-tailed Student-t blocks.

The folder layout mirrors the C++ suite: one subfolder per demo, plus a single
consolidated DGP file.

```txt
trex_da/
  ├── README.md
  ├── trex_da_dgps.R                   <- all scenario DGPs (one file)
  ├── demo_trex_da_01_ar1/
  │   ├── demo_trex_da_01_ar1.R
  │   ├── README.md
  │   └── simulation_results/
  ├── demo_trex_da_02_equi/
  │   └── ...
  └── demo_trex_da_08_groups/
```

Dependencies: `TRexSelectorNeo`, `plotly` (correlation heatmaps),
`parallel` (MC loops). Each demo sources the consolidated
[trex_da_dgps.R](trex_da_dgps.R) (all scenario generators in one file, mirroring
the C++ `dgp_generators.hpp`) and the shared
[../support_generators.R](../support_generators.R) /
[../simulation_utils.R](../simulation_utils.R).

---

## Demos

| File | Description | cpp counterpart |
|---|---|---|
| [demo_trex_da_01_ar1](demo_trex_da_01_ar1/) | AR(1) Toeplitz DGP: single run; SNR sweep (Random vs CappedSpread support); rho sweep; 2D gap x rho sweep exploring the kappa boundary | `demo_trex_da_01_mc_sim_ar1` |
| [demo_trex_da_02_equi](demo_trex_da_02_equi/) | Full-matrix equicorrelated DGP (BT factor-model DGP with `rho_within == rho_between`) | `demo_trex_da_02_mc_sim_equi_and_bt` |
| [demo_trex_da_03_nn](demo_trex_da_03_nn/) | Nearest-Neighbours / MA(kappa) banded DGP: single run; SNR, rho, and kappa sweeps; 2D kappa x rho sweep | `demo_trex_da_03_mc_sim_nn` |
| [demo_trex_da_03b_nn_ar](demo_trex_da_03b_nn_ar/) | TREX+DA+NN applied to AR(1) data (method-mismatch test): single run; SNR, rho, and SNR x rho sweeps | `demo_trex_da_03b_mc_sim_nn_ar` |
| [demo_trex_da_04_bt_ar1_block](demo_trex_da_04_bt_ar1_block/) | Pure AR(1)-block DGP (p = M*Q): SNR, rho, Q, and M sweeps, each with an outer linkage loop | `demo_trex_da_04_mc_sim_bt_ar1_block` |
| [demo_trex_da_05_bt_ar1_block_white](demo_trex_da_05_bt_ar1_block_white/) | AR(1)-block + white-noise DGP (fixed p_total, p_white = p_total − M*Q): SNR, rho, Q, and M sweeps with linkage loop | `demo_trex_da_05_mc_sim_bt_ar1_block_sweeps` |
| [demo_trex_da_06_bt_ht_block](demo_trex_da_06_bt_ht_block/) | Heavy-tailed Student-t(3) AR(1)-block data, Gaussian vs Student-t noise scenarios: SNR, rho, Q, M, tFDR, and linkage sweeps | `demo_trex_da_06_mc_sim_bt_ht_block_sweeps` |
| [demo_trex_da_07_bt_ht_block_white](demo_trex_da_07_bt_ht_block_white/) | Heavy-tailed blocks + heavy-tailed white augmentation, Gaussian vs heavy-tailed noise: SNR, rho, Q, M, and tFDR sweeps | `demo_trex_da_07_mc_sim_bt_ht_block_white` |
| [demo_trex_da_08_groups](demo_trex_da_08_groups/) | Prior-groups method: 3-level hierarchy (groups of 10 / 50 / 250) over a multi-level latent-factor DGP with Toeplitz leaf blocks, SNR sweep. The method is selected by supplying a non-empty `prior_groups` to `trex_da_control()` | `demo_trex_da_08_mc_sim_groups` |

Notes:

- Demos 06 and 07 redirect their console output to a log file via `sink()`
  (comment out the `sink()` lines to print to the console instead).
- The prior-groups method (demo 08) is driven by passing `prior_groups` and
  `rho_grid_labels` to `trex_da_control()`; the C++ core routes on
  `prior_groups` being non-empty, so `da_method` is then irrelevant. (The R
  binding gained these arguments in the 2026-07 upstream fix, alongside
  `cv_seed` / `cv_n_folds` / `cv_n_lambda` for GVS.)

---

## DGPs

All scenario generators live in the single [trex_da_dgps.R](trex_da_dgps.R)
(mirroring the C++ `dgp_generators.hpp`):

- `dgp_ar1_snr` — AR(1) Toeplitz design, SNR-calibrated noise.
- `dgp_bt_snr` / `dgp_ar1_block_snr` / `dgp_ar1_block_white_snr` — binary-tree /
  hierarchical-block and AR(1)-block covariance (with/without white augmentation).
- `dgp_nn_snr` / `dgp_nn_block_snr` — banded MA(kappa) covariance via causal
  convolution of a shared innovation sheet.
- `dgp_ht_block_snr` / `dgp_ht_block_white_snr` — heavy-tailed Student-t AR(1)
  blocks via Gaussian scale mixture.
- `dgp_groups_toeplitz_leaf` — multi-level latent-factor design with Toeplitz
  leaf blocks (prior-groups demo).

All DGPs calibrate the noise variance to a target SNR = Var(signal)/Var(noise)
and take the true support as 1-based R indices. Earlier demo-less generators
(the `dgp_arp_snr` AR(p) family and the `dgp_hapgen_snr` LD-block DGP) were
removed in the consolidation.

---

## Recorded results

Manually-recorded result summaries live in each demo's `simulation_results/`
folder, e.g.
[demo_trex_da_01_ar1/simulation_results/data/Results_trex_da_01_ar1.md](demo_trex_da_01_ar1/simulation_results/data/Results_trex_da_01_ar1.md).

---

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_01_ar1/demo_trex_da_01_ar1.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_01_ar1/demo_trex_da_01_ar1.R 8    # 8 worker cores
```

Scripts resolve their own directory, so they run from any working directory.
Several demos gate individual parts behind `if (FALSE)` / run flags — open the
file and enable the parts you want. Each demo folder has its own `README.md`.

---

**Last updated**: 2026-07-08
