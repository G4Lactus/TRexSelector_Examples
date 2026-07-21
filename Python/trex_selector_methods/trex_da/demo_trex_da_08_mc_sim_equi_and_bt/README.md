# Demo 08: DA-TRex+EQUI on Equicorrelated (EQUI) Data — 2D SNR x rho Sweep (Python)

Monte-Carlo simulations for **DA-TRex+EQUI** vs. the **classical (no-DA) T-Rex
selector** on equicorrelated data. One 2D sweep: for each SNR column in
`{0.5, 1, 2, 5}`, a full rho sweep over
`{0, 0.025, 0.05, 0.1, 0.2, 0.4, 0.7, 0.9}` — dense near 0 to locate the DA
suppression cliff. Common configuration: n=300, p=1000, s=10, amplitude 3.0,
tFDR=0.2, K=20 random experiments, MC=200 per grid point, `Random` support with
a full redraw per trial; solvers TLARS / TAFS / TOMP, each with its base
(no-DA) T-Rex comparison row — 6 method rows per grid point.

Python port of the C++ demo
[cpp/trex_selector_methods/trex_da/demo_trex_da_08_mc_sim_equi_and_bt](../../../../cpp/trex_selector_methods/trex_da/demo_trex_da_08_mc_sim_equi_and_bt/README.md)
— same DGP, control parameters, sweep axes, console summaries, and result-file
naming (the folder keeps the cpp name "equi_and_bt": earlier revisions also
carried BT sweeps on block-equicorrelated data, dropped as redundant since the
BT method is exercised on better-suited block designs in Demos 02-05).

The greedy solvers use *exchangeable tie-breaking* (`exch_tie_alpha = 0.25` for
TAFS/TOMP, `0` for TLARS); TAFS additionally runs with `rho_afs = 0.3`.

> **Purpose**: equicorrelation is the stress case of the `trex_da` framework —
> every column loads on one shared factor, so there is no "distance" to
> exploit the way AR(1) has, and individual variable selection is structurally
> infeasible for any method. The 2D layout shows that (1) the DA suppression
> cliff sits at very small rho, (2) it is SNR-independent, and (3) at rho=0
> everything works. DA-EQUI fails *safe* (selects nothing); base T-Rex fails
> *unsafe* (FDR far above target).

## Data-generating process

`dgp_equi_snr` generates rows through a single latent factor,
`X[i,j] = sqrt(rho) * f_i + sqrt(1 - rho) * eta[i,j]`, giving the compound
symmetry covariance `(1 - rho) I + rho 11^T`. As rho -> 1 the p-1 minor
eigenvalues collapse and the regression becomes unidentifiable; unlike AR(1),
correlation is uniform across all variable pairs, so the shared factor puts
every variable inside every other variable's DA correction set. Support:
Random, redrawn per trial; SNR-calibrated noise.

## Sweep design

- **Outer loop (SNR columns)**: `SNR in {0.5, 1, 2, 5}`, one output file pair
  per column, disjoint seed blocks per column (statistically independent).
- **Inner sweep (rho)**: `{0, 0.025, 0.05, 0.1, 0.2, 0.4, 0.7, 0.9}` —
  front-loaded at small rho where the suppression cliff is expected.

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_08_mc_sim_equi_and_bt/demo_trex_da_08_mc_sim_equi_and_bt.py
```

Monte Carlo trials run in parallel via `ProcessPoolExecutor` (up to
`min(6, os.cpu_count())` workers). Reduce `NUM_MC` to preview.

## Output files

Written to `simulation_results/` (`.txt` pretty-printed, `.csv` tidy long
format), one pair per SNR column:

- `da_trex_mc_da_equi_rho_snr0p5.txt` / `.csv`
- `da_trex_mc_da_equi_rho_snr1.txt` / `.csv`
- `da_trex_mc_da_equi_rho_snr2.txt` / `.csv`
- `da_trex_mc_da_equi_rho_snr5.txt` / `.csv`

## What to expect

See the C++ counterpart README for the committed reference results: the DA
suppression cliff sits at the first nonzero grid point (rho=0.025) in every
SNR column — signal strength never moves it; base T-Rex TLARS fails unsafe
from the same grid point (FDR 0.36-0.40, climbing to 0.92-0.97 at rho=0.9);
at SNR=5 base TOMP is the only usable selector for moderate equicorrelation.
Rare DA-TOMP lone-survivor false positives at high rho in the zero-power
regime are a confirmed extreme-order-statistic mechanism, not a bug.

---

**Last updated**: 2026-07-21
