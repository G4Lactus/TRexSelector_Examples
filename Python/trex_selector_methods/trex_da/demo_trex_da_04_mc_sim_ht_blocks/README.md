# Demo 04: DA-TRex+BT on Heavy-Tailed Block Data (Python)

Monte-Carlo simulations for **DA-TRex+BT** — the Binary-Tree Dependency-Aware
T-Rex selector (`DAMethod.BT`) — on a **heavy-tailed** block-diagonal Toeplitz
design: the predictors always follow a multivariate Student-t(nu=3) law, and
the response noise is either Gaussian (**s1_Gauss**) or t(nu=3) (**s2_Heavy**).
Five ceteris-paribus sweeps (SNR, rho, Q, M, target FDR) plus a dedicated
linkage sweep, each run under both noise scenarios and the three HAC linkage
methods. Common setup: n=150, amplitude 1.0, rho=0.8, nu=3, tFDR=0.2, K=20
random experiments, MC=200 per grid point; solvers TLARS / TAFS / TOMP.

Python port of the C++ demo
[cpp/trex_selector_methods/trex_da/demo_trex_da_04_mc_sim_ht_blocks](../../../../cpp/trex_selector_methods/trex_da/demo_trex_da_04_mc_sim_ht_blocks/README.md)
— same DGP, control parameters, sweep axes, console summaries, and result-file
naming.

The greedy solvers use *exchangeable tie-breaking* (`exch_tie_alpha = 0.25` for
TAFS/TOMP, `0` for TLARS) in Parts 1-4; the cpp tFDR and linkage sections set
solver_type / rho_afs / stagnation but not `exch_tie_alpha`, and Parts 5-6
mirror that (the greedy solvers run there with `exch_tie_alpha = 0`). TAFS
runs with `rho_afs = 0.3` throughout.

## Data-generating process

`dgp_block_toeplitz_hvt_snr` builds M independent AR(1)-Toeplitz blocks of
size Q (`p = M*Q`, no white padding) and converts each block to a
**multivariate Student-t(nu)** via a Gaussian scale mixture with an independent
chi-squared mixing variable per (row, block) — extreme events are block-local.
The `OnePerBlock` support places one active per block (s = M, amplitude 1.0).
The two scenarios differ only in the noise law: `s1_Gauss` is Gaussian,
`s2_Heavy` is variance-matched t(nu) — the predictors are heavy-tailed in
*both* scenarios. Base point: M=5, Q=5 (p=25, s=5), rho=0.8, nu=3, SNR=2.0,
seed 2026.

## What it runs

Parts are toggled by the `RUN_PART_*` flags (all `True`). Every non-linkage
part loops scenario x linkage.

- **Part 1**: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- **Part 2**: rho sweep `{0.0, 0.1, ..., 0.9}`.
- **Part 3**: Q sweep `{5, 10, ..., 50}`; `p = M*Q` grows to 250 > n.
- **Part 4**: M sweep `{1, 3, 5, 10, 15, 20, 30}`; `s = M`, `p = 5M`.
- **Part 5**: tFDR sweep `{0.05, 0.10, ..., 0.50}` (realized-FDR calibration).
- **Part 6**: linkage sweep, encoded `{1=Single, 2=Complete, 3=Average}`.

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_04_mc_sim_ht_blocks/demo_trex_da_04_mc_sim_ht_blocks.py
```

Monte Carlo trials run in parallel via `ProcessPoolExecutor` (up to
`min(6, os.cpu_count())` workers). Reduce `NUM_MC` to preview.

## Output files

Written to `simulation_results/` (64 files — the largest output set of any
demo here: 32 scenario stems, one `.txt`+`.csv` pair each):

- `da_trex_mc_da_ht_blocks_{snr,rho,Q,M,tFDR}_{s1_Gauss,s2_Heavy}_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ht_blocks_linkage_{s1_Gauss,s2_Heavy}.txt` / `.csv`

## What to expect

See the C++ counterpart README for the committed reference results: TOMP is the
clear winner (high TPR, lowest FDR); TLARS grazes and then breaches the
tFDR=0.2 target from SNR~0.5 on and sits above target across the whole Q sweep
— heavy-tailed predictors, not dimension growth, loosen its FDR control;
Gaussian vs. heavy-tailed *noise* is nearly indistinguishable because the
heavy-tailed predictors dominate; the M sweep collapses power for M >~ 15-20;
linkage choice barely matters at the well-correlated base point.

---

**Last updated**: 2026-07-21
