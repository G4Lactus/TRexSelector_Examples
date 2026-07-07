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

Dependencies: `TRexSelectorNeo`, `plotly` (correlation heatmaps),
`parallel` (MC loops). The demos source the shared
[../support_generators.R](../support_generators.R) and
[../simulation_utils.R](../simulation_utils.R) plus the local `dgp_*.R` files.

---

## Demos

| File | Description | cpp counterpart |
|---|---|---|
| [demo_trex_da_01_ar1.R](demo_trex_da_01_ar1.R) | AR(1) Toeplitz DGP: single run; SNR sweep (Random vs CappedSpread support); rho sweep; 2D gap x rho sweep exploring the kappa boundary | `demo_trex_da_01_mc_sim_ar1` |
| [demo_trex_da_02_equi.R](demo_trex_da_02_equi.R) | Full-matrix equicorrelated DGP (BT factor-model DGP with `rho_within == rho_between`) | `demo_trex_da_02_mc_sim_equi_and_bt` |
| [demo_trex_da_03_bt_equi_blocks.R](demo_trex_da_03_bt_equi_blocks.R) | Binary-Tree (BT) hierarchical-block model with equicorrelated blocks: single run; SNR sweep; linkage sweep over single/complete/average | `demo_trex_da_02_mc_sim_equi_and_bt` |
| [demo_trex_da_04_nn_01_ma.R](demo_trex_da_04_nn_01_ma.R) | Nearest-Neighbours / MA(kappa) banded DGP: single run; SNR, rho, and kappa sweeps; 2D kappa x rho sweep | `demo_trex_da_03_mc_sim_nn` |
| [demo_trex_da_04_nn_02_ar.R](demo_trex_da_04_nn_02_ar.R) | TREX+DA+NN applied to AR(1) data (method-mismatch test): single run; SNR, rho, and SNR x rho sweeps | `demo_trex_da_03b_mc_sim_nn_ar` |
| [demo_trex_da_05_bt_ar1_block_sweeps.R](demo_trex_da_05_bt_ar1_block_sweeps.R) | Pure AR(1)-block DGP (p = M*Q): SNR, rho, Q, and M sweeps, each with an outer linkage loop | `demo_trex_da_04_mc_sim_bt_ar1_block` |
| [demo_trex_da_06_bt_ar1_plus_white_block_sweeps.R](demo_trex_da_06_bt_ar1_plus_white_block_sweeps.R) | AR(1)-block + white-noise DGP (fixed p_total, p_white = p_total − M*Q): SNR, rho, Q, and M sweeps with linkage loop | `demo_trex_da_05_mc_sim_bt_ar1_block_sweeps` |
| [demo_trex_da_07_bt_heavy_tailed_sweeps.R](demo_trex_da_07_bt_heavy_tailed_sweeps.R) | Heavy-tailed Student-t(3) AR(1)-block data, Gaussian vs Student-t noise scenarios: SNR, rho, Q, M sweeps | `demo_trex_da_06_mc_sim_bt_ht_block_sweeps` |
| [demo_trex_da_08_bt_heavy_tailed_plus_white_block_sweeps.R](demo_trex_da_08_bt_heavy_tailed_plus_white_block_sweeps.R) | Heavy-tailed blocks + identity-covariance Student-t augmentation, Gaussian vs heavy-tailed noise: SNR, rho, Q, M, and tFDR sweeps | `demo_trex_da_07_mc_sim_bt_ht_block_white` |

Notes:

- Demos 07 and 08 redirect their console output to a log file via `sink()`
  (comment out the `sink()` lines to print to the console instead).
- The C++ suite additionally contains `demo_trex_da_08_mc_sim_groups`
  (multi-level nested group structures), which has no R counterpart here.

### Known limitation

`trex_da_control(da_method = ...)` in the installed **TRexSelectorNeo** R
binding accepts only `"AR1"`, `"EQUI"`, and `"BT"`. The C++ core supports
`DAMethod::NN`, but the R binding does not expose `"NN"` yet. The two
nearest-neighbour demos (`demo_trex_da_04_nn_01_ma.R`,
`demo_trex_da_04_nn_02_ar.R`) therefore carry a guard that prints a `[SKIP]`
notice and exits cleanly under `Rscript` until the binding adds `"NN"`.

---

## DGP helper files

| File | Description |
|---|---|
| [dgp_ar1_snr.R](dgp_ar1_snr.R) | `dgp_ar1_snr()` — AR(1) Toeplitz design generated column-by-column, SNR-calibrated noise |
| [dgp_arp_snr.R](dgp_arp_snr.R) | AR(p) generalization with Yule-Walker-calibrated unit-variance columns; coefficient families `make_ar_phi_geometric()` / `make_ar_phi_pacf()`; AR(1) case reproduces `dgp_ar1_snr()` exactly (not sourced by any current demo) |
| [dgp_bt_snr.R](dgp_bt_snr.R) | Binary-tree / hierarchical-block DGPs: `dgp_bt_snr()` (factor-model block covariance), `dgp_ar1_block_snr()`, `dgp_ar1_block_white_snr()` |
| [dgp_ht_snr.R](dgp_ht_snr.R) | Heavy-tailed Student-t AR(1)-block DGPs via Gaussian scale mixture: `dgp_ht_block_snr()`, `dgp_ht_block_white_snr()`; Gaussian or Student-t noise |
| [dgp_nn_snr.R](dgp_nn_snr.R) | `dgp_nn_snr()` — banded MA(kappa) covariance via causal convolution of a shared innovation sheet |
| [dgp_hapgen_snr.R](dgp_hapgen_snr.R) | Analytic LD-block correlation structure mirroring `hapgen_like.py` (6 blocks, p = 150 SNPs; not sourced by any current demo) |

All DGPs calibrate the noise variance to a target
SNR = Var(signal) / Var(noise) and take the true support as 1-based R indices.

---

## Recorded results

Completed simulation runs are documented as markdown reports in
[simulation_results/](simulation_results/):

- [Results_trex_da_01_ar1.md](simulation_results/Results_trex_da_01_ar1.md)
- [Results_trex_da_02_equi.md](simulation_results/Results_trex_da_02_equi.md)
- [Results_trex_da_03_bt_equi_blocks.md](simulation_results/Results_trex_da_03_bt_equi_blocks.md)
- [Results_trex_da_04_nn.md](simulation_results/Results_trex_da_04_nn.md)
- [Results_trex_da_05_bt_ar1_block_sweeps.md](simulation_results/Results_trex_da_05_bt_ar1_block_sweeps.md)
- [Results_trex_da_06_bt_ar1_pwhiteblock_sweeps.md](simulation_results/Results_trex_da_06_bt_ar1_pwhiteblock_sweeps.md)
- [Results_trex_da_07_bt_ht_sweeps.md](simulation_results/Results_trex_da_07_bt_ht_sweeps.md)
- [Results_trex_da_08_bt_ht_pwhiteblock_sweeps.md](simulation_results/Results_trex_da_08_bt_ht_pwhiteblock_sweeps.md)

---

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_01_ar1.R
```

Scripts resolve their own directory, so they run from any working directory.
Several demos gate individual parts behind `if (FALSE)` / run flags — open the
file and enable the parts you want.

---

**Last updated**: 2026-07-08
