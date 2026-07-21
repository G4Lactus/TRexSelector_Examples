# Demo 01: Screen-TRex on an i.i.d. Gaussian Design (Python)

## Purpose

Establishes the **baseline behaviour of Screen-TRex** on the easiest possible
design — independent Gaussian predictors with no correlation structure. Two
thresholding rules are compared over a signal-to-noise sweep: the **Ordinary**
majority-vote rule (`Phi_j > 0.5`) and the **Bootstrap-CI** rule. Everything the
later demos add (correlation, dependency-aware corrections, biobank routing,
solver backends) is measured against the numbers established here. In-memory
variant — see Demo 02 for the memory-mapped counterpart. Python port of
[`cpp/trex_selector_methods/trex_screening/demo_trex_scr_01_mc_sim_screen_trex`](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_01_mc_sim_screen_trex/README.md)
(see that README for the statistical model and the annotated result discussion).

## DGP / Data

`dgp_iid_snr` — i.i.d. standard-Gaussian design with SNR-calibrated noise,
`sigma = sqrt(||X beta||^2 / (n * SNR))`.

- `n = 300`, `p = 1000` (high-dimensional, `p > n`), `s = 10`.
- The active support is drawn **uniformly at random per Monte Carlo trial**
  (RNG offset `trial_seed + 500000`, isolated from the design/noise draws), with
  coefficients `beta_j = 1`.
- SNR sweep: `{0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.

## Methods compared

Two thresholding rules, both `ScreenTRexMethod.TREX`:

- **Screen-TRex Ordinary: TLARS** — `bootstrap=False`.
- **Screen-TRex Bootstrap-CI: TLARS** — `bootstrap=True` (`R_boot = 1000`
  replicates, `ci_grid_step = 0.001`).

Shared solver control: `solver = "TLARS"`, `K = 20`, `use_memory_mapping = False`.
Each trial builds `tsn.TRexScreeningSelector(X, y, screen_control=…,
trex_control=…, seed=-1)` and reads `sel.selected_indices` (0-based) and
`res.estimated_FDR`.

## What it computes

`demo_screen_trex_monte_carlo()` sweeps the two methods x the SNR grid, calling
`run_mc_screen` per cell (`base_seed = 24 + snr_idx * 1000`, mirroring the cpp
source). It records mean realized FDR, TPR, and the procedure's own estimated
FDR. `num_MC = 100` (committed default; the cpp demo uses 200 — downscaled for a
practical Python runtime, override with `SCR_NUM_MC`).

## Imports & parallelism

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `trex_scr_common` module; the package imports as
`trex_selector_neo`. The MC loop is parallelized with
`concurrent.futures.ProcessPoolExecutor` (`SCR_NUM_WORKERS`, default 6).

## Output

Written to `simulation_results/data/`:

- `scr_screen_trex_snr_n300_p1000.txt` / `.csv` — FDR, TPR, and estimated FDR
  per method and SNR level (tidy long CSV: `method, metric, SNR, value`).

Until the demo is run, `simulation_results/` holds only a `.gitkeep`.

## Running

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_01_mc_sim_screen_trex/demo_trex_scr_01_mc_sim_screen_trex.py
```

## Interpretation

- Screening is **not FDR-controlling at low SNR**: at `SNR <= 0.2` almost every
  selection is a false one (TPR near zero), a property of thresholding a voting
  statistic without a target level.
- Both rules become reliable from `SNR >= 0.5`; the Ordinary rule recovers all
  ten active variables at `SNR = 5`.
- **Bootstrap-CI is the conservative rule**: lower FDR, clearly less power, and
  slower per run due to the bootstrap replicates.
- The internal FDR estimate is **conservative here** (sits above the realized
  FDR once signal is present) — Demo 03 shows it failing dramatically under
  correlation.

**Last updated**: 2026-07-21
