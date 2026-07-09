# Demo 05: T-Rex+GVS on Heavy-Tailed (Student-t(3)) Blocks (Python)

## Purpose

Stress-test T-Rex+GVS on non-Gaussian data. This demo reuses the same four-block
active/trap geometry as Demo 03, but replaces every Gaussian component of the
data-generating process — latent factors, within-block perturbations, and response
noise — with scaled Student-t(3) variables, testing robustness to heavy-tailed
departures from the Gaussian assumption.

Python port of `cpp/.../demo_trex_gvs_05_mc_sim_t3_blocks`. The R counterpart is
`R/trex_selector_methods/trex_gvs/demo_trex_gvs_05/`. Column indices are **0-based**
here (Python convention); the R counterpart is 1-based.

## Data-generating process

`dgp_t3_equi_blocks_snr` (from `trex_gvs_dgps.py`) — four equicorrelated blocks of
sizes `{20, 50, 80, 65}`: three active (`beta = 3`, `s = 150`) plus one inactive
heavy-tailed trap, shuffled into `p = 500` with white-noise gaps. Columns follow
`X_j = sqrt(rho) * Z_g + sqrt(1 - rho) * E`, where the factors `Z`, the
perturbations `E`, and the response noise are drawn from `t(df) / sqrt(df/(df-2))`
with `df = 3` (rescaled to unit variance). Within-block equicorrelation `rho` is
passed directly to the DGP. Dimensions: `n = 200`, `p = 500`. SNR-calibrated noise.

## Methods compared

- **EN** — `gvs_type = "EN"`, `en_solver = "TENET"`
- **EN+AUG** — `gvs_type = "EN"`, `en_solver = "TENET_AUG"`
- **IEN** — `gvs_type = "IEN"`, `en_solver = "TENET"`

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "single"`,
`df = 3`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — Single-run demo: heavy-tail check (`max |X|` should exceed the
  Gaussian ~3), block layout (shuffled order + block-order ID sequence), a
  within-block correlation check on block 0 vs `rho`, and one EN / EN+AUG / IEN run.
- **Part 2** — SNR sweep at fixed `rho = 0.75`; `SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`.
- **Part 3** — rho sweep at fixed `SNR = 2.0`;
  `rho in {0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99}`.
- **Part 4** — 2-D SNR x rho grid (5 SNR x 6 rho, rho axis `{0.30, 0.50, 0.70, 0.85, 0.90, 0.99}`).

Unlike the R Part 1, the Python Part 1 does **not** render a Plotly correlation
heatmap; it prints the heavy-tail / block-layout / within-block-correlation
diagnostics above and then runs the three methods.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_05/demo_trex_gvs_05.py       # Part 1 only (default)
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_05/demo_trex_gvs_05.py 8     # 8 worker cores
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

- EN and EN+AUG should remain reasonably FDP-controlled with TPP rising toward the
  mid-0.9s as SNR grows, though slightly below the Gaussian equicorrelated case
  (Demo 03) at matched SNR — heavy tails inject occasional large outliers that blur
  group boundaries.
- Unlike the Gaussian equicorrelated demos, IEN is expected to be much less
  conservative here: FDP tightening below target as SNR grows while TPP climbs into
  the 0.7–0.8 range — a more competitive pattern.
- Key question over the SNR sweep: does heavy-tailed noise mainly cost power (lower
  TPP), or also compromise FDP control (drifting above 0.1)?

**Last updated**: 2026-07-08
