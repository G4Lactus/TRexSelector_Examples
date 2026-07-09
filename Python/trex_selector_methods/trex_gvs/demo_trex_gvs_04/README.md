# Demo 04: T-Rex+GVS on Negative-Correlation Traps (Python)

## Purpose

Test whether T-Rex+GVS can recover an **active group with negative within-group
correlation** (sign-flipped halves) while **excluding two spatially separated
inactive trap groups** that share the same sign-flipped pattern. This stresses the
grouping logic's ability to distinguish signal from confounding correlation
geometry, not just correlation strength.

Python port of `cpp/.../demo_trex_gvs_04_mc_sim_neg_traps`. The R counterpart is
`R/trex_selector_methods/trex_gvs/demo_trex_gvs_04/`. Column indices are **0-based**
here (Python convention); the R counterpart is 1-based.

## Data-generating process

`dgp_neg_traps_snr` (from `trex_gvs_dgps.py`) — each correlated group loads half its
members `+Z` and half `-Z`, inducing negative within-group correlation. All ranges
are 0-based:

- **Group 1** (cols 0–99): ACTIVE, `beta = +3` for cols 0–49 and `beta = -3` for
  cols 50–99, `s = 100`.
- **Trap 1** (cols 100–199): inactive, same sign-flipped structure.
- **Noise** (cols 200–299): white noise.
- **Trap 2** (cols 300–359): second, smaller inactive sign-flipped trap.
- **Noise** (cols 360–499): white noise.

Within-group correlation follows `rho = 1/(1 + sd_x^2)`. Dimensions: `n = 200`,
`p = 500`. Default `sd_x = sqrt(0.01)` (`rho ~ 0.99`). SNR-calibrated noise.

## Methods compared

- **EN** — `gvs_type = "EN"`, `en_solver = "TENET"`
- **EN+AUG** — `gvs_type = "EN"`, `en_solver = "TENET_AUG"`
- **IEN** — `gvs_type = "IEN"`, `en_solver = "TENET"`

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "single"`,
`num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — Single-run demo: negative-correlation structural checks
  (`Cor(X[:,0], X[:,50])` in the active group and `Cor(X[:,300], X[:,330])` in
  Trap 2, both expected `~ -rho`), one EN / EN+AUG / IEN run, and a false-positive
  autopsy classifying each FP as leaked from Trap 1 (cols 100–199), Trap 2
  (cols 300–359), or white noise.
- **Part 2** — SNR sweep at fixed `sd_x = sqrt(0.01)` (`rho ~ 0.99`);
  `SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`.
- **Part 3** — rho sweep at fixed `SNR = 2.0` (`sd_x` derived from rho);
  `rho in {0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99}`.
- **Part 4** — 2-D SNR x rho grid (5 SNR x 6 rho).

Unlike the R Part 1, the Python Part 1 does **not** render a Plotly correlation
heatmap; it prints the negative-correlation diagnostics above and then runs the
three methods.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_04/demo_trex_gvs_04.py       # Part 1 only (default)
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_04/demo_trex_gvs_04.py 8     # 8 worker cores
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

- EN and EN+AUG should track each other closely with FDP near the `tFDR = 0.1`
  target and TPP reaching near-complete recovery at moderate SNR — the sign-flipped
  active group is still identified as one group via the (unsigned) correlation-based
  HAC step, and both traps should be excluded.
- IEN is expected to stay conservative here (low FDP, slowly rising TPP), similar to
  the uniformly-positive-correlation scenarios; the sign flips alone do not change
  IEN behavior much.
- The two spatially separated traps are a stronger test of false-group exclusion
  than a single trap — watch the FP autopsy and whether FDP creeps up with SNR.

**Last updated**: 2026-07-08
