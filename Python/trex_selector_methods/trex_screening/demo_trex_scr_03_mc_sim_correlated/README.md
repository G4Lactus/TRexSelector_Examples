# Demo 03: Screen-TRex under Correlated Designs (Python)

## Purpose

Measures what **predictor correlation does to Screen-TRex**, and whether the
**dependency-aware (DA)** variants repair the damage. Three correlation
structures are put in front of the screener — AR(1), equicorrelated, and
block-equicorrelated — and in each case the plain Ordinary and Bootstrap-CI
rules are compared against their DA counterparts, which pre-group correlated
variables before the votes are counted. Python port of
[`cpp/trex_selector_methods/trex_screening/demo_trex_scr_03_mc_sim_correlated`](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_03_mc_sim_correlated/README.md)
(see that README for the correlation models and the detailed per-part results).

## DGP / Data

`n = 300`, `p = 1000`, `p1 = 10` evenly-spaced active entries, columns
standardised to unit variance so `rho` is the realized column correlation. Three
generators: `dgp_ar1` (`corr = rho^|j-j'|`), `dgp_equi` (one shared latent
factor), `dgp_block_equi` (one latent factor per block, 5 contiguous blocks).
The design and noise are drawn afresh in every trial.

- SNR grid: `{0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- rho grid: `{0.1, 0.2, …, 0.9}`.

## Methods compared

Every part compares four methods — the two thresholding rules, plain and
dependency-aware: **Ordinary: TLARS**, **Bootstrap-CI: TLARS**,
**Ordinary+DA-\*: TLARS**, **Bootstrap-CI+DA-\*: TLARS**, where the DA variant
(`TREX_DA_AR1` / `TREX_DA_EQUI` / `TREX_DA_BLOCK_EQUI`) is always matched to the
design that generated the data. DA grouping uses `rho_thr_DA = 0.02`; the
block-equi part passes `n_blocks = 5` to the screen control.

## The five parts (matching the cpp source)

| Part | Design | Fixed | Sweep | DA variant | Output stem |
| :--- | :----- | :---- | :---- | :--------- | :---------- |
| 1 | AR(1) | `rho = 0.5`, `beta = 1` | SNR | DA-AR1 | `scr_ar1_snr_n300_p1000_rho0.50` |
| 2 | Equicorrelated | `rho = 0.4`, `beta = 3` | SNR | DA-EQUI | `scr_equi_snr_n300_p1000_rho0.40` |
| 3 | AR(1) | `SNR = 1`, `beta = 1` | rho | DA-AR1 | `scr_ar1_rho_n300_p1000_snr1.00` |
| 4 | Equicorrelated | `SNR = 1`, `beta = 3` | rho | DA-EQUI | `scr_equi_rho_n300_p1000_snr1.00` |
| 5 | Block-equi (5 blocks) | `SNR = 1`, `beta = 3` | rho | DA-BLOCK-EQUI | `scr_block_equi_rho_n300_p1000_blocks5_snr1.00` |

Each part calls `run_mc_screen` per (method, grid point), `base_seed = 42 +
idx * 1000`. `num_MC = 50` (committed default; cpp uses 200 — downscaled across
five heavy MC parts, override with `SCR_NUM_MC`). The MC loop is parallelized
with `ProcessPoolExecutor` (`SCR_NUM_WORKERS`, default 6).

## Output

Each part writes `.txt` + tidy `.csv` (`method, metric, SNR|rho, value`) to
`simulation_results/data/` under the stems above.

## Running

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_03_mc_sim_correlated/demo_trex_scr_03_mc_sim_correlated.py
```

## Interpretation

- **AR(1) correlation is survivable**; DA-AR1 lowers the error rate at a modest
  power cost and buys roughly one step of the rho grid.
- **Equicorrelation breaks screening outright** — even at moderate `rho = 0.4`
  the FDR sits at 0.8–0.96 — and no DA variant fully rescues the Ordinary rule.
- **Block structure sits between the two extremes**, closer to the bad one; DA
  helps meaningfully only at low correlation.
- The central caution: under equicorrelation the **estimated FDR is wildly
  optimistic exactly where the method fails**, off by a factor of twenty or
  more — the opposite of the safe-side behaviour of Demo 01.

**Last updated**: 2026-07-21
