# Demo 02: T-Rex+DA+EQUI on Equicorrelated Data (Python)

## Purpose

Cover the equicorrelation (compound-symmetry) dependency structure with the `EQUI`
DA method. A single-run illustration is paired with Monte Carlo sweeps over SNR and
rho.

Python port of `R/trex_selector_methods/trex_da/demo_trex_da_02_equi` and the equi
parts of `cpp/.../demo_trex_da_02_mc_sim_equi_and_bt`.

## Data-generating process

Reuses the two-level factor-model DGP `dgp_bt_snr` with `rho_within == rho_between`,
which collapses to pure equicorrelation: `X[i,j] = sqrt(rho) * f_i + sqrt(1-rho) *
eta[i,j]`. Single-run config: `n=150, p=500, s=5, n_blocks=5, rho=0.25, snr=2.0,
amplitude=3.0, tFDR=0.2, K=20, seed=2026`, Random support. Monte Carlo config:
`n=300, p=1000, s=10, amplitude=3.0, n_blocks=10, tFDR=0.2, K=20, num_MC=200,
seed=2026`, Random support redrawn per trial.

## DA method

`build_da_selector("EQUI", ...)` constructs a `trex_selector_neo.TRexDASelector`
with a `TRexDAControlParameter` (`method = DAMethod.EQUI`, no extra DA arguments) and
a `TRexControlParameter` (`K = 20`, `solver_type = SolverTypeForTRex.TLARS`).

## What it runs

Parts are toggled by the `RUN_PART_*` flags. Default: Part 1 only
(`RUN_PART_1 = True`; Parts 2-3 `False`).

- Part 1: single-run demo on equicorrelated data (rho=0.25, SNR=2.0); prints selection and TPP/FDP.
- Part 2: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` at rho=0.25.
- Part 3: 1D pure-equi rho sweep `{0.0, 0.1, ..., 0.9}` at SNR=2.0 (`rho_within == rho_between`).

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_da/demo_trex_da_02_equi/demo_trex_da_02_equi.py
```

The file inserts both its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` and `da_sim_common` modules import directly. Enable the
Monte Carlo parts by flipping the `RUN_PART_*` flags. Output is console-only; the
folder carries its own `simulation_results/` directory for optionally recorded
summaries. Monte Carlo parts run in parallel via Python `multiprocessing` using the
`spawn` start method (up to `min(6, os.cpu_count())` workers).

## Interpretation

- Because equicorrelation ties every column to the same shared factor, the EQUI correction behaves qualitatively like the AR(1) case in Demo 01 (reduced FDR vs base T-Rex at some TPR cost), but more uniformly across variables, since equicorrelation has no notion of column distance.

**Last updated**: 2026-07-08
