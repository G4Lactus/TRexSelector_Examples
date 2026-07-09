# Demo 04: T-Rex+GVS on Negative-Correlation Traps

## Purpose

Test whether T-Rex+GVS can recover an **active group with negative within-group correlation** (sign-flipped halves) while **excluding two spatially separated inactive trap groups** that share the same sign-flipped pattern. This stresses the grouping logic's ability to distinguish signal from confounding correlation geometry, not just correlation strength.

## Data-generating process

`dgp_neg_traps_snr` — each correlated group loads half its members `+Z` and half `-Z`, inducing negative within-group correlation:

- **Group 1** (cols 1–100): ACTIVE, `beta = +3 / -3`, `s = 100`.
- **Trap 1** (cols 101–200): inactive, same sign-flipped structure.
- **Noise** (cols 201–300): white noise.
- **Trap 2** (cols 301–360): second, smaller inactive sign-flipped trap.
- **Noise** (cols 361–500): white noise.

Within-group correlation follows `rho = 1/(1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`. Default `sd_x = sqrt(0.01)` (`rho ~ 0.99`).

## Methods compared

- **EN** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **EN+AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **IEN** — `GVS_type = "IEN"`, `en_solver = "TENET"`

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — Single-run demo: negative-correlation checks (expected `~ -rho` within active group and trap 2), one EN / EN+AUG / IEN run, and a false-positive autopsy classifying each FP as leaked from Trap 1, Trap 2, or white noise.
- **Part 2** — SNR sweep at fixed `sd_x = sqrt(0.01)` (`rho ~ 0.99`); `SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`.
- **Part 3** — rho sweep at fixed `SNR = 2.0`; `rho in {0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99}`.
- **Part 4** — 2-D SNR x rho grid (5 SNR x 6 rho).

## Running

```bash
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_04/demo_trex_gvs_04.R      # all parts
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_04/demo_trex_gvs_04.R 8    # 8 worker cores
```

Output is console-only (no result files are written); the `run_part_*` flags at the top toggle individual parts; the script sources the shared `../trex_gvs_dgps.R` and `../../simulation_utils.R`.

## Interpretation

- EN and EN+AUG should track each other closely with FDP near the `tFDR = 0.1` target and TPP reaching near-complete recovery at moderate SNR — the sign-flipped active group is still identified as one group via the (unsigned) correlation-based HAC step, and both traps should be excluded.
- IEN is expected to stay conservative here (low FDP, slowly rising TPP), similar to the uniformly-positive-correlation scenarios; the sign flips alone do not change IEN behavior much.
- The two spatially separated traps are a stronger test of false-group exclusion than a single trap — watch the FP autopsy and whether FDP creeps up with SNR.

**Last updated**: 2026-07-08
