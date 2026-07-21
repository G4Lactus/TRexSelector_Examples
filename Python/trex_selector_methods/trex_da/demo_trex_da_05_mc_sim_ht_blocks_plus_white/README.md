# Demo 05: DA-TRex+BT on Heavy-Tailed Block + White-Noise Data (Python)

Monte-Carlo simulations for **DA-TRex+BT** — the Binary-Tree Dependency-Aware
T-Rex selector (`DAMethod.BT`) — on the "combined effects" design: heavy-tailed
(multivariate t(nu)) AR(1)-Toeplitz blocks embedded in a larger design padded
with i.i.d. Student-t white-noise columns. It merges Demo 02's white-noise
dilution with Demo 04's heavy-tailed robustness study, sweeping SNR, rho, Q, M,
and the target level tFDR — each across two noise scenarios (Gaussian vs.
heavy-tailed response noise) and the three HAC linkage methods.
Common setup: n=150, p_total=500, amplitude 1.0, nu=3, tFDR=0.2 (except the
tFDR sweep), K=20 random experiments, MC=200 per grid point; solvers
TLARS / TAFS / TOMP.

Python port of the C++ demo
[cpp/trex_selector_methods/trex_da/demo_trex_da_05_mc_sim_ht_blocks_plus_white](../../../../cpp/trex_selector_methods/trex_da/demo_trex_da_05_mc_sim_ht_blocks_plus_white/README.md)
— same DGP, control parameters, sweep axes, console summaries, and result-file
naming.

The greedy solvers use *exchangeable tie-breaking* (`exch_tie_alpha = 0.25` for
TAFS/TOMP, `0` for TLARS) in Parts 1-4; the cpp tFDR section sets
solver_type / rho_afs / stagnation but not `exch_tie_alpha`, and Part 5 mirrors
that. TAFS runs with `rho_afs = 0.3` throughout.

## Data-generating process

`dgp_ht_block_white_snr` concatenates M independent **heavy-tailed** blocks of
size Q with a white-noise block. Each block starts from the Q x Q AR(1)
Toeplitz correlation and is transformed into a multivariate t(nu) vector via a
per-(row, block) chi-squared scale mixture — when a block draws a mixing
variable near zero, all Q variables of that block explode simultaneously
(block-local tail dependence). The remaining `p_white = p_total - M*Q` columns
are i.i.d. t(nu) white noise. The `OnePerBlock` support places one active per
heavy-tailed block (s = M, amplitude 1.0). Noise scenarios: `s1_Gauss`
(Gaussian) vs. `s2_Heavy` (variance-matched t(nu)); SNR-calibrated in both.
Base point: M=5, Q=5 (p_ar=25, p_white=475, s=5), rho=0.8, nu=3, SNR=2.0,
seed 2026.

## What it runs

Parts are toggled by the `RUN_PART_*` flags (all `True`). Every part loops
scenario x linkage.

- **Part 1**: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- **Part 2**: rho sweep `{0.0, 0.1, ..., 0.9}`.
- **Part 3**: Q sweep `{5, 10, ..., 50}`; `p_ar = M*Q`, `p_white = 500 - p_ar`.
- **Part 4**: M sweep `{1, 3, 5, 10, 15, 20, 30}`; `p_ar`, `p_white`, `s = M`
  all vary.
- **Part 5**: tFDR sweep `{0.05, 0.10, ..., 0.50}`.

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_05_mc_sim_ht_blocks_plus_white/demo_trex_da_05_mc_sim_ht_blocks_plus_white.py
```

Monte Carlo trials run in parallel via `ProcessPoolExecutor` (up to
`min(6, os.cpu_count())` workers). Reduce `NUM_MC` to preview.

## Output files

Written to `simulation_results/` (60 files = 30 scenario stems, one
`.txt`+`.csv` pair each):

- `da_trex_mc_da_ht_blocks_plus_white_{snr,rho,Q,M,tFDR}_{s1_Gauss,s2_Heavy}_{Single,Complete,Average}.txt` / `.csv`

## What to expect

See the C++ counterpart README for the committed reference results: the
white-noise padding again *helps* FDR control (well below Demo 04's levels in
the Gaussian scenario), while heavy-tailed *response* noise pushes cells to
(and marginally over) the target at high SNR / high rho; the M-sweep power
collapse hits earlier and harder than in the Gaussian demos; realized FDR
tracks at or below the target across the tFDR sweep. TOMP is the strongest
solver, TAFS unusable in the single-signal M=1 corner.

---

**Last updated**: 2026-07-21
