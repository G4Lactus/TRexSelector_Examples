# Demo 02: DA-TRex+BT on Block AR(1) + White-Noise Data (Python)

Monte-Carlo simulations for **DA-TRex+BT** — the Binary-Tree Dependency-Aware
T-Rex selector (`DAMethod.BT`) — on a block-diagonal AR(1) design diluted by a
large block of i.i.d. white-noise columns, sweeping SNR, rho, block size Q, and
number of blocks M, each across the three HAC linkage methods
(Single / Complete / Average).
Common setup: n=300, p_total=1000, amplitude 1.0, tFDR=0.1, K=20 random
experiments, MC=200 per grid point; solvers TLARS / TAFS / TOMP.

Python port of the C++ demo
[cpp/trex_selector_methods/trex_da/demo_trex_da_02_mc_sim_ar1_blocks_plus_white](../../../../cpp/trex_selector_methods/trex_da/demo_trex_da_02_mc_sim_ar1_blocks_plus_white/README.md)
— same DGP, control parameters, sweep axes, console summaries, and result-file
naming.

The greedy solvers use *exchangeable tie-breaking* (`exch_tie_alpha = 0.25` for
TAFS/TOMP, `0` for TLARS); TAFS additionally runs with `rho_afs = 0.3`.

## Data-generating process

`dgp_ar1_block_white_snr` concatenates M statistically independent AR(1)
blocks of size Q (within-block Toeplitz correlation `rho^|j-k|`, generated via
the AR(1) recursion; `p_ar = M*Q`) with `p_white = p_total - p_ar` i.i.d.
N(0,1) columns — the white block is pure null padding that makes the problem
high-dimensional (`p_total > n`). The `OnePerBlock` support places exactly one
active variable per AR block (s = M, all actives in the AR part, amplitude
1.0); noise is SNR-calibrated. Base point: M=5, Q=5 (p_ar=25, p_white=975,
s=5), rho=0.7, SNR=2.0, seed 2026.

The **BT group design** clusters variables by HAC on the dissimilarity
`1 - |corr(x_j, x_j')|`; the dendrogram cut defines disjoint groups for the
nearest-partner occurrence deflation. The **linkage method** (Single /
Complete / Average) sets the between-cluster distance — whether that choice
matters on a clean block design is this demo's core question, so every sweep
runs once per linkage.

## What it runs

Parts are toggled by the `RUN_PART_*` flags (all `True`).

- **Part 1**: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- **Part 2**: rho sweep `{0.0, 0.1, ..., 0.9}`.
- **Part 3**: Q sweep `{5, 10, ..., 50}`; `p_ar = M*Q` varies,
  `p_white = 1000 - p_ar`; support recomputed per grid point.
- **Part 4**: M sweep `{1, 3, 5, 10, 15, 20, 30}`; `p_ar`, `p_white`, and
  `s = M` all vary.

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_02_mc_sim_ar1_blocks_plus_white/demo_trex_da_02_mc_sim_ar1_blocks_plus_white.py
```

Monte Carlo trials run in parallel via `ProcessPoolExecutor` (up to
`min(6, os.cpu_count())` workers). Reduce `NUM_MC` to preview.

## Output files

Written to `simulation_results/` (24 files = 12 scenario stems, one
`.txt`+`.csv` pair each):

- `da_trex_mc_da_ar1_blocks_plus_white_snr_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_blocks_plus_white_rho_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_blocks_plus_white_Q_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_blocks_plus_white_M_{Single,Complete,Average}.txt` / `.csv`

## What to expect

See the C++ counterpart README for the committed reference results: FDR is very
tightly controlled (max ~0.044 on the SNR sweep — the 975 uncorrelated null
columns dilute the design); Single linkage loses power at *low* rho and large
Q, while Complete/Average are robust; the M sweep shows a graceful, FDR-safe
power collapse once M >~ 20 as the fixed signal budget spreads over more active
variables. Compare against Demo 03 (same block core, no white padding) to
isolate the effect of the null dilution.

---

**Last updated**: 2026-07-21
