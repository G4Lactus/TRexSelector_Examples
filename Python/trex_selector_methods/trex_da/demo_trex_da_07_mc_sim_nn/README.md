# Demo 07: DA-TRex+NN on Banded Data + AR(1) Method-Mismatch Study (Python)

Monte-Carlo simulations for **DA-TRex+NN** — the Nearest-Neighbor
Dependency-Aware T-Rex selector (`DAMethod.NN`) — on banded MA(kappa) data
(the correctly-specified case) plus a **method-mismatch stress test** applying
the same NN correction to AR(1) data.
Common setup: n=300, p=1000, s=10, amplitude 3.0, tFDR=0.2, K=20 random
experiments, MC=200 per grid point; solvers TLARS / TAFS / TOMP (Part 1
additionally runs the classical no-DA T-Rex baseline).

Python port of the C++ demo
[cpp/trex_selector_methods/trex_da/demo_trex_da_07_mc_sim_nn](../../../../cpp/trex_selector_methods/trex_da/demo_trex_da_07_mc_sim_nn/README.md)
— same DGPs, control parameters, sweep axes, console summaries, and
result-file naming. The AR(1) mismatch study was its own demo (03b) in the
legacy Python numbering; it is Part 3 here.

The greedy solvers use *exchangeable tie-breaking* (`exch_tie_alpha = 0.25` for
TAFS/TOMP, `0` for TLARS) in Part 1; the cpp 2D loops of Parts 2-3 set
solver_type / rho_afs / stagnation but not `exch_tie_alpha`, and the Python
port mirrors that. TAFS runs with `rho_afs = 0.3` throughout.

## Data-generating processes

**NN / MA(kappa) model (Parts 1-2, `dgp_nn_snr`).** Each predictor is a causal
MA(kappa) convolution of a shared innovation sheet with geometric coefficients
`theta_l ~ rho^l`, normalised to unit variance. The covariance is banded
Toeplitz: correlation decays roughly geometrically inside the band and cuts
off exactly at bandwidth kappa — the hallmark difference to AR(1), whose
correlation never truly vanishes.

**AR(1) model (Part 3, `dgp_ar1_snr` — same DGP as Demo 01).** Geometrically
decaying but unbanded correlation — deliberately mismatched to the NN
correction's finite-range group design.

Supports: `CappedSpread` (deterministic even spacing, max_gap=20) and `Random`
(redrawn per trial); SNR-calibrated noise.

## What it runs

Parts are toggled by the `RUN_PART_*` flags (all `True`).

- **Part 1**: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` on NN data (kappa=3,
  rho=0.7), CappedSpread and Random support, with the base (no-DA) T-Rex
  comparison rows.
- **Part 2**: 2D kappa x rho sweep on NN data;
  `kappa in {1, 2, 3, 5, 7, 10, 15}` x `rho in {0.1, 0.3, 0.5, 0.7, 0.8, 0.9}`
  at SNR=2.0; CappedSpread and Random sub-sections per solver.
- **Part 3**: 2D SNR x rho sweep on AR(1) data with DA+NN (method mismatch);
  `SNR in {0.1, ..., 5}` x `rho in {0.1, ..., 0.9}`.

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_07_mc_sim_nn/demo_trex_da_07_mc_sim_nn.py
```

Monte Carlo trials run in parallel via `ProcessPoolExecutor` (up to
`min(6, os.cpu_count())` workers). Reduce `NUM_MC` to preview.

## Output files

Written to `simulation_results/` (`.txt` pretty-printed, `.csv` tidy long
format):

- `da_trex_mc_da_nn_snr_capped.txt` / `.csv`,
  `da_trex_mc_da_nn_snr_random.txt` / `.csv` (Part 1)
- `da_trex_mc_da_nn_kappa_rho.txt` / `.csv` (Part 2, two-support 2D schema)
- `da_trex_mc_da_nn_ar_snr_rho.txt` / `.csv` (Part 3, two-support 2D schema)

## What to expect

See the C++ counterpart README for the committed reference results: on NN data
all three DA solvers stay below tFDR=0.2 at every SNR while base T-Rex TLARS
violates from SNR=0.5 on; in the kappa x rho grid, FDR violations appear only
in the top-right corner (rho >= 0.8 and kappa >= 7). On the AR(1) mismatch,
control holds up to rho ~0.7 and breaks down sharply at rho = 0.8-0.9 (up to
~0.65 for TLARS) — the price of the wrong dependency model, negligible at
moderate correlation and large in the strong-correlation regime (compare
Demo 01, where the correctly-specified AR1 correction stays controlled).

---

**Last updated**: 2026-07-21
