# Demo 03: DA-TRex+BT on Block-Diagonal AR(1) Data (Python)

Monte-Carlo simulations for **DA-TRex+BT** — the Binary-Tree Dependency-Aware
T-Rex selector (`DAMethod.BT`) — on a clean block-diagonal AR(1) design (no
white-noise padding), sweeping SNR, rho, block size Q, and number of blocks M,
each across the three HAC linkage methods (Single / Complete / Average).
Common setup: n=150, amplitude 1.0, tFDR=0.2, K=20 random experiments, MC=200
per grid point; solvers TLARS / TAFS / TOMP.

Python port of the C++ demo
[cpp/trex_selector_methods/trex_da/demo_trex_da_03_mc_sim_ar1_blocks](../../../../cpp/trex_selector_methods/trex_da/demo_trex_da_03_mc_sim_ar1_blocks/README.md)
— same DGP, control parameters, sweep axes, console summaries, and result-file
naming.

The greedy solvers use *exchangeable tie-breaking* (`exch_tie_alpha = 0.25` for
TAFS/TOMP, `0` for TLARS); TAFS additionally runs with `rho_afs = 0.3`.

## Data-generating process

`dgp_ar1_block_snr` builds M statistically independent AR(1) blocks of size Q —
and nothing else: unlike Demo 02 there is no white-noise padding, so every
column belongs to a correlated block and `p = M*Q`. Within-block correlation is
Toeplitz `rho^|j-k|` (AR(1) recursion per block); the `OnePerBlock` support
places one active variable per block (s = M, amplitude 1.0); noise is
SNR-calibrated. Base point: M=5, Q=5 (p=25, s=5), rho=0.7, SNR=2.0, seed 2026.
At the base point the problem is low-dimensional (p=25 << n=150); the Q sweep
pushes it past p > n (p=250 at Q=50) and the M sweep reaches p = n = 150 at
M=30.

The BT group design and the linkage question are as in
[Demo 02](../demo_trex_da_02_mc_sim_ar1_blocks_plus_white/README.md).

## What it runs

Parts are toggled by the `RUN_PART_*` flags (all `True`).

- **Part 1**: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- **Part 2**: rho sweep `{0.0, 0.1, ..., 0.9}`.
- **Part 3**: Q sweep `{5, 10, ..., 50}`; `p = M*Q` varies; s = M = 5 fixed.
- **Part 4**: M sweep `{1, 3, 5, 10, 15, 20, 30}`; `p = M*Q` and `s = M` vary.

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_03_mc_sim_ar1_blocks/demo_trex_da_03_mc_sim_ar1_blocks.py
```

Monte Carlo trials run in parallel via `ProcessPoolExecutor` (up to
`min(6, os.cpu_count())` workers). Reduce `NUM_MC` to preview.

## Output files

Written to `simulation_results/` (24 files = 12 scenario stems, one
`.txt`+`.csv` pair each):

- `da_trex_mc_da_ar1_blocks_snr_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_blocks_rho_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_blocks_Q_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_blocks_M_{Single,Complete,Average}.txt` / `.csv`

## What to expect

See the C++ counterpart README for the committed reference results: FDR stays
under target on the SNR sweep but visibly higher than the white-noise-diluted
Demo 02 (every false-discovery candidate is a correlated shadow here); on the
rho sweep TLARS breaks the target at rho=0.9 under every linkage (~0.25) while
TOMP stays controlled; Single linkage again loses power at low rho and large
Q; the M sweep collapses power for M >~ 20 while FDR stays below target (the
selector fails safe). The M=1 corner is degenerate — TAFS selects nothing at
all with a single signal.

---

**Last updated**: 2026-07-21
