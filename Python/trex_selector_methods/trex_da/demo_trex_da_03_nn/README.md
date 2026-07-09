# Demo 03: T-Rex+DA+NN on Banded (Nearest-Neighbor) Data (Python)

## Purpose

Study the `NN` DA method, which targets a banded MA(kappa) covariance where only
the kappa nearest-neighboring columns are correlated (`Sigma[j,k] != 0` only for
`|j-k| <= kappa`) — a sharp, finite-range dependency in contrast to AR(1)'s
geometric decay across all distances.

Python port of `R/trex_selector_methods/trex_da/demo_trex_da_04_nn_01_ma` and
`cpp/.../demo_trex_da_03_mc_sim_nn`.

## Data-generating process

`dgp_nn_snr` builds a banded MA(kappa) design where columns j and k are correlated
only when `|j-k| <= kappa`; the MA coefficients decay geometrically in rho.
Single-run config: `n=150, p=500, s=5, kappa=1, rho=0.7, snr=2.0, amplitude=3.0,
tFDR=0.2, K=20, seed=2026`. Monte Carlo config: `n=300, p=1000, s=10, amplitude=3.0,
kappa=3, rho=0.7, snr=2.0, tFDR=0.2, K=20, num_MC=200, seed=2026, max_gap=20`. Sweeps
use CappedSpread(max_gap=20) support.

## DA method

`build_da_selector("NN", ...)` constructs a `trex_selector_neo.TRexDASelector` with
a `TRexDAControlParameter` (`method = DAMethod.NN`, no extra DA arguments) and a
`TRexControlParameter` (`K = 20`, `solver_type = SolverTypeForTRex.TLARS`).

## What it runs

Parts are toggled by the `RUN_PART_*` flags. Default: Parts 1 and 6 active
(`RUN_PART_1 = True`, `RUN_PART_6 = True`; Parts 2-5 `False`), matching the R demo.

- Part 1: single-run demo (kappa=1, rho=0.7, SNR=2.0); prints selection and TPP/FDP.
- Part 2: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` (kappa=3, rho=0.7).
- Part 3: rho sweep `{0.1, ..., 0.9}` (kappa=3, SNR=2.0).
- Part 4: kappa sweep `{1, 2, 3, 5, 7, 10, 15}` (rho=0.7, SNR=2.0).
- Part 5: 2D kappa x rho sweep; kappa `{1,2,3,5,7,10,15}` x rho `{0.1,0.3,0.5,0.7,0.8,0.9}` at SNR=2.0.
- Part 6: 2D SNR x rho sweep; SNR `{0.1,0.5,1.0,2.0,5.0}` x rho `{0.5,0.6,0.7,0.75,0.8,0.85,0.9}` at kappa=3.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_da/demo_trex_da_03_nn/demo_trex_da_03_nn.py
```

The file inserts both its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` and `da_sim_common` modules import directly. Toggle
parts via the `RUN_PART_*` flags. Output is console-only; the folder carries its own
`simulation_results/` directory for optionally recorded summaries. Monte Carlo parts
run in parallel via Python `multiprocessing` using the `spawn` start method (up to
`min(6, os.cpu_count())` workers).

## Interpretation

- Because NN correlation has a hard cutoff at distance kappa, the correction is most effective (best FDR control, least TPR cost) when active variables are spaced farther apart than kappa.
- As kappa grows the effective correlation neighborhood widens and NN needs progressively wider correction windows, with a corresponding TPR cost. Compare against Demo 03b, which applies the same NN correction to AR(1) data (mismatched).

**Last updated**: 2026-07-08
