# Demo 03: Screen-TRex on Correlated Designs

## Purpose

The centerpiece correlation-robustness demo: compares **Ordinary** screening against **dependency-aware (DA)** screening variants — `TREX_DA_AR1`, `TREX_DA_EQUI`, `TREX_DA_BLOCK_EQUI` — on AR(1), equicorrelated, and block-equicorrelated designs, where naive Ordinary screening is expected to lose power because correlated variables split voting evidence.

> **Status**: `simulation_results/` is currently empty — this demo has not yet been run in this checkout.

---

## Data Generation Parameters

All nine parts (four single-run + five Monte Carlo) use $n=300$, $p=1000$, $p_1=10$ active variables (evenly spaced via `make_beta()`), with correlation structure depending on the part:

- **AR(1)**: $X_{i,j} = \rho\, X_{i,j-1} + \sqrt{1-\rho^2}\, w_{ij}$ (columns standardized to unit variance afterward).
- **Equicorrelated**: shared-latent-factor construction (see `dgp_equi`).
- **Block-equicorrelated**: `n_blocks` equal-size correlated blocks (see `dgp_block_equi`).

> **Note**: the source file's top-of-file doc-comment and the inline comments inside `main()` number these nine parts inconsistently with each other. The tables below identify each part by its DGP/method content rather than by a "Demo N" number to avoid repeating that ambiguity.

### Single-run parts (local seed 1, selector seed 42)
| Method | DGP | $\rho$ | SNR | $\beta$ |
|---|---|---|---|---|
| DA-AR1 | AR(1) | 0.5 | 1.0 | 1.0 |
| DA-EQUI | Equicorrelated | 0.4 | 1.0 | 1.0 |
| Ordinary (no DA, baseline) | AR(1) | 0.7 | 5.0 | 3.0 |
| DA-BLOCK-EQUI | Block-equicorrelated (`n_blocks=5`) | 0.5 | 1.0 | 1.0 |

### Monte Carlo parts (MC = 500 per grid point)
| Sweep | DGP | Fixed parameter | Swept values |
|---|---|---|---|
| SNR sweep | AR(1) | $\rho=0.5$, $\beta=1.0$ | $\mathrm{SNR}\in\{0.01,0.1,0.2,0.5,0.6,1.0,2.0,5.0\}$ |
| SNR sweep | Equicorrelated | $\rho=0.4$, $\beta=3.0$ | same SNR grid |
| $\rho$ sweep | AR(1) | $\mathrm{SNR}=1.0$, $\beta=1.0$ | $\rho \in \{0.1,\dots,0.9\}$ |
| $\rho$ sweep | Equicorrelated | $\mathrm{SNR}=1.0$, $\beta=3.0$ | $\rho \in \{0.1,\dots,0.9\}$ |
| $\rho$ sweep | Block-equicorrelated (`n_blocks=5`) | $\mathrm{SNR}=1.0$, $\beta=3.0$ | $\rho \in \{0.1,\dots,0.9\}$ |

Each MC sweep part compares **Ordinary**, **Bootstrap-CI**, and the matching **DA** variant.

---

## Control Parameters

```
K = 20
solver = TLARS
rho_thr_DA = 0.02
n_blocks = 5   # block-equicorrelated part only
MC = 500 per grid point
```

---

## Output Files (expected)

Written to `simulation_results/` once run (exact stems from the source code):

- `d03_screen_trex_da_ar1_mc_n300_p1000_rho0.50.txt` / `.csv` — SNR sweep, AR(1).
- `d03_screen_trex_da_equi_mc_n300_p1000_rho0.40.txt` / `.csv` — SNR sweep, equicorrelated.
- `d03_screen_trex_da_ar1_rho_sweep_n300_p1000_snr1.00.txt` / `.csv` — $\rho$ sweep, AR(1).
- `d03_screen_trex_da_equi_rho_sweep_n300_p1000_snr1.00.txt` / `.csv` — $\rho$ sweep, equicorrelated.
- `d03_screen_trex_da_block_equi_rho_sweep_n300_p1000_blocks5_snr1.00.txt` / `.csv` — $\rho$ sweep, block-equicorrelated.

The four single-run parts (DA-AR1, DA-EQUI, Ordinary-on-AR1 baseline, DA-BLOCK-EQUI) print diagnostics to the console only and do not write files.

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_screening/demo_trex_scr_03_screen_trex_correlated/demo_trex_scr_03_screen_trex_correlated
```

---

## Interpretation

- The single-run **Ordinary-on-AR(1)** part is explicitly labeled a baseline "no DA correction" comparison — it is meant to be read against the DA-AR1 single-run part on the same correlation family to see the qualitative effect of dependency-aware screening at a glance.
- In the SNR-sweep and $\rho$-sweep MC parts, expect **Ordinary** screening to degrade (lower TPR, and/or estimated FDR drifting away from realized FDR) as $\rho$ increases, while the matching **DA** variant should hold up better across the $\rho$ grid — this is the central hypothesis this demo is designed to test.
- **Bootstrap-CI** (without DA correction) is a useful middle comparison: check whether bootstrap-based thresholding alone recovers some of the power lost by Ordinary screening under correlation, or whether the DA pre-grouping step is doing the real work.
- The block-equicorrelated part additionally tests whether `n_blocks` (fixed at 5 here) is well matched to the DGP's own block count — a mismatch here is conceptually analogous to the `corr_max` sensitivity studied for HAC clustering in the `trex_gvs` demos.

---

**Last updated**: 2026-07-08
