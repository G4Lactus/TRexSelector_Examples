# Demo 06: T-Rex+GVS on AR(1) Blocks (Python)

## Purpose

Evaluate T-Rex+GVS when within-block correlation follows an **AR(1) time-series
structure** instead of equicorrelation. This uses the same four-block active/trap
geometry as Demo 03, but replaces the shared-latent-factor model with an
AR(1)/Toeplitz covariance within each block, testing robustness to realistic,
decaying correlation patterns.

Python port of `cpp/.../demo_trex_gvs_06_mc_sim_ar1_blocks`. The R counterpart is
`R/trex_selector_methods/trex_gvs/demo_trex_gvs_06/`. Column indices are **0-based**
here (Python convention); the R counterpart is 1-based.

## Data-generating process

`dgp_ar1_blocks_snr` (from `trex_gvs_dgps.py`) — four blocks of sizes
`{20, 50, 80, 65}`: three active (`beta = 3`, `s = 150`) plus one inactive AR(1)
trap, shuffled into `p = 500` with i.i.d. `N(0,1)` gaps. Within-block correlation is
AR(1)/Toeplitz: `Cor(X_j, X_k) = rho^|j-k|`, so correlation decays with column
distance and a block no longer behaves like one tight cluster. The AR(1) parameter
`rho` is passed directly to the DGP. Dimensions: `n = 200`, `p = 500`.
SNR-calibrated noise.

## Methods compared

- **EN** — `gvs_type = "EN"`, `en_solver = "TENET"`
- **EN+AUG** — `gvs_type = "EN"`, `en_solver = "TENET_AUG"`
- **IEN** — `gvs_type = "IEN"`, `en_solver = "TENET"`

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "single"`,
`num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — Single-run demo: block layout (shuffled order + block-order ID
  sequence), an AR(1) decay check on block 0 (empirical lag-1/2/3 correlations vs
  theoretical `rho^k`), and one EN / EN+AUG / IEN run.
- **Part 2** — SNR sweep at fixed `rho = 0.85`; `SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`.
- **Part 3** — rho sweep at fixed `SNR = 2.0`;
  `rho in {0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99}`.
- **Part 4** — 2-D SNR x rho grid (5 SNR x 6 rho, rho axis `{0.30, 0.50, 0.70, 0.85, 0.90, 0.99}`).

Unlike the R Part 1, the Python Part 1 does **not** render a Plotly correlation
heatmap; it prints the AR(1)-decay and block-layout diagnostics above and then runs
the three methods.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_06/demo_trex_gvs_06.py       # Part 1 only (default)
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_06/demo_trex_gvs_06.py 8     # 8 worker cores
```

Only **Part 1 runs by default** (`RUN_PART_1 = True`); the heavy Monte Carlo parts
(200 MC per grid point) are **opt-in** via the `RUN_PART_2/3/4` flags at the top of
the file. This differs from the R suite, where all parts default on. An optional
first CLI argument sets the worker-core count (otherwise `min(6, os.cpu_count())`).

The file inserts both its own directory and the parent suite directory onto
`sys.path`, so the shared `trex_gvs_dgps` / `gvs_sim_common` modules resolve and are
re-importable by the `spawn`-launched Monte Carlo workers — so the demo must be run
as a file, not piped via `python -`. Output is console-only (no result files are
written).

## Interpretation

- FDP should stay low (often well below `tFDR = 0.1`) across the SNR grid for all
  three methods, but TPP should recover much more slowly with SNR than in the
  equicorrelated Hastie/Mixed-Blocks demos — AR(1) decay means only nearby columns
  are strongly correlated, weakening group-level evidence aggregation.
- IEN is expected to sit between its near-zero TPP on the equicorrelated Gaussian
  demos and the EN-based methods, showing real (if smaller) power gains as SNR
  increases.
- The rho sweep is the most informative view: watch how quickly recovery degrades as
  `rho` decreases and the AR(1) decay steepens.

**Last updated**: 2026-07-08
