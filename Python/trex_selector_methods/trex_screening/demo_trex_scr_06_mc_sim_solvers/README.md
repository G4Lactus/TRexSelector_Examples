# Demo 06: Solver Backends for Screen-TRex — TLARS, TAFS, TOMP (Python)

## Purpose

Asks whether the **choice of T-Rex solver backend** changes screening
performance, or whether it is merely a speed knob. The same i.i.d. Gaussian
design as [Demo 01](../demo_trex_scr_01_mc_sim_screen_trex/README.md) is used, so
the numbers are directly comparable. Three backends are swept over SNR — **TLARS**
(terminating LARS), **TAFS-0.3** (terminating adaptive forward selection,
`rho_afs = 0.3`), and **TOMP** (terminating orthogonal matching pursuit) — first
under the Ordinary rule alone, then crossed with Bootstrap-CI. Python port of
[`cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_mc_sim_solvers`](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_mc_sim_solvers/README.md).

## DGP / Data

Same DGP as Demo 01: `dgp_iid_snr`, `n = 300`, `p = 1000`, `s = 10`, random
support per trial (`beta_j = 1`), SNR sweep
`{0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}` — any difference in the curves is
attributable to the solver backend, not the data.

## Methods compared

Three backends driving `ScreenTRexMethod.TREX`, each built via a flat control
that (mirroring the cpp `make_trex_ctrl_for()`) enables the T-loop stagnation
stop (`tloop_stagnation_stop = True`, `tloop_max_stagnant_steps = 5`) for every
backend:

- **TLARS** — `solver = "TLARS"`.
- **TAFS-0.3** — `solver = "TAFS"`, `solver_params.rho_afs = 0.3`.
- **TOMP** — `solver = "TOMP"`.

## The two parts (matching the cpp source)

- **Part 1 — solvers under the Ordinary rule.** Three series
  (`Screen-TRex Ordinary: TLARS / TAFS-0.3 / TOMP`) over the SNR grid. Stem
  `scr_solvers_snr_n300_p1000_s10`.
- **Part 2 — solvers crossed with both rules.** Six series
  (`Ordinary` and `Bootstrap-CI` prefixes x the three backends). Stem
  `scr_solver_method_snr_n300_p1000_s10`.

Each part calls `run_mc_screen` per (variant, SNR), `base_seed = 24 +
snr_idx * 1000`. `num_MC = 50` (committed default; cpp uses 200 — downscaled,
override with `SCR_NUM_MC`). Trials run in a `ProcessPoolExecutor`
(`SCR_NUM_WORKERS`, default 6).

## Output

Written to `simulation_results/data/` as `.txt` + tidy `.csv`
(`method, metric, SNR, value`) under the two stems above.

## Running

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_06_mc_sim_solvers/demo_trex_scr_06_mc_sim_solvers.py
```

## Interpretation

- The solver is **not a pure speed knob**: at mid SNR the greedier backends
  (TAFS-0.3, TOMP) are markedly more powerful than TLARS at comparable FDR; all
  three converge at high SNR.
- **Under Bootstrap-CI the backend matters far more than under the Ordinary
  rule** — TLARS is the most conservative, TOMP the most powerful.
- Practical takeaway: TLARS is the safe conservative default; TOMP and TAFS-0.3
  buy substantial mid-SNR power at a comparable error rate.

**Last updated**: 2026-07-21
