# Demo 03: Screen-TRex on Correlated Designs

## Purpose

The centerpiece correlation-robustness demo: compares **Ordinary** screening against the **dependency-aware (DA)** variants — `trex_screen_control(trex_method = "TREX_DA_AR1" / "TREX_DA_EQUI" / "TREX_DA_BLOCK_EQUI")` — on AR(1), equicorrelated, and block-equicorrelated designs. Ordinary screening thresholds each `Phi_j > 0.5` independently, so redundant correlated variables split the voting evidence and power collapses; the DA variants pre-group strongly correlated variables before voting and additionally report an estimated rho via `sel$estimated_correlation`.

**Important (R binding)**: `trex_screen_control()` exposes only `use_bootstrap_CI`, `R_boot`, `ci_grid_step`, and `trex_method`. Unlike the C++ `make_screen_control()`, it does **not** expose `rho_thr_DA` or `n_blocks` — the DA variants use the library's internal defaults for those, so the block DGP's `n_blocks = 5` is not forwarded to the screener.

## Data-generating processes

All nine parts use `n = 300`, `p = 1000`, `p1 = 10` active variables (evenly spaced via `make_beta()`), SNR-calibrated noise. From the shared `../trex_scr_sim_utils.R`:

- `dgp_ar1` — AR(1) (banded) correlation.
- `dgp_equi` — equicorrelation (single shared latent factor).
- `dgp_block_equi` — block-equicorrelation (`n_blocks = 5` equal-size blocks, per-block factors).

Grids: `SNR in {0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`, `rho in {0.1, 0.2, ..., 0.9}`.

## Methods compared

Each MC sweep runs four `TRexScreeningSelector` configurations: **Ordinary** (`trex_method = "TREX"`), **Bootstrap-CI** (`"TREX"` + `use_bootstrap_CI = TRUE`), and the matching **DA** variant with and without bootstrap (e.g. Ord+DA-AR1 / Boot+DA-AR1). Shared controls via `make_scr_trex_ctrl()`: `solver = "TLARS"`, `K = 20`; selector `seed = 42`.

`num_MC = 10` per grid point (the C++ demo uses 500; reduced so all five MC parts finish in about two minutes). Results are statistically comparable to the C++ reference, **not** bit-identical — the two use different RNG streams. Selected indices and support sets are 1-based in R (C++ is 0-based).

## Parts

Single runs (Parts 1–4, DGP seed 1, console-only):

- **Part 1** — Ordinary on AR(1) (`rho = 0.7`, `SNR = 5`, `beta = 3`): baseline, no DA correction.
- **Part 2** — DA-AR1 on AR(1) (`rho = 0.5`, `SNR = 1`, `beta = 1`).
- **Part 3** — DA-EQUI on equicorrelated (`rho = 0.4`, `SNR = 1`, `beta = 1`).
- **Part 4** — DA-BLOCK-EQUI on block-equicorrelated (`rho = 0.5`, 5 blocks, `SNR = 1`, `beta = 1`).

Monte Carlo sweeps (Parts 5–9, four methods each):

- **Part 5** — SNR sweep on AR(1) (`rho = 0.5`, `beta = 1`) → `d03_da_ar1_mc_n300_p1000_rho0.50`.
- **Part 6** — SNR sweep on equicorrelated (`rho = 0.4`, `beta = 3`) → `d03_da_equi_mc_n300_p1000_rho0.40`.
- **Part 7** — rho sweep on AR(1) (`SNR = 1`, `beta = 1`) → `d03_da_ar1_rho_sweep_n300_p1000_snr1.00`.
- **Part 8** — rho sweep on equicorrelated (`SNR = 1`, `beta = 3`) → `d03_da_equi_rho_sweep_n300_p1000_snr1.00`.
- **Part 9** — rho sweep on block-equicorrelated (`SNR = 1`, `beta = 3`, 5 blocks) → `d03_da_block_equi_rho_sweep_n300_p1000_snr1.00`.

## Running

```bash
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_03_screen_correlated/demo_trex_scr_03_screen_correlated.R
```

The script sources the shared `../trex_scr_sim_utils.R`; the `run_single_parts` / `run_mc_parts` flags at the top toggle the two halves. Each MC sweep writes a human-readable `.txt` and a tidy `.csv` (FDR / TPR / Estimated FDR per method and grid point) into the local `simulation_results/` folder; Parts 1–4 print to the console only.

## Interpretation

- Read Part 1 (Ordinary on AR(1), no DA correction) against Part 2 (DA-AR1 on the same correlation family) to see the qualitative effect of dependency-aware screening at a glance.
- In the recorded AR(1) rho sweep (Part 7), Ordinary screening's realized FDR climbs with correlation (about 0.05–0.12 for `rho <= 0.5`, rising to about 0.53 at `rho = 0.9`) while its Estimated FDR stays flat near 0.1 — the self-reported FDR becomes anti-conservative under strong correlation. Ord+DA-AR1 keeps realized FDR lower across most of the grid at comparable TPR.
- Bootstrap-CI (without DA) is the middle comparison: in the recorded runs it lowers FDR but at a substantial TPR cost, suggesting the DA pre-grouping, not bootstrap thresholding, does the real work; Boot+DA combinations are very conservative (near-zero FDR, low TPR).
- The equicorrelated rho sweep (Part 8) is the hardest recorded setting: all four methods show high realized FDR at moderate-to-high rho with Estimated FDR far below it — plausibly because the R binding cannot pass `rho_thr_DA` / `n_blocks`, so DA grouping runs on internal defaults. Treat these curves as a binding-limitation illustration rather than the method's C++-tuned behavior.
- With only 10 MC repetitions per grid point, individual cells are noisy (see e.g. the Bootstrap dip at `rho = 0.6` in Part 7); trends across the grid are meaningful, single values are not.

**Last updated**: 2026-07-08
