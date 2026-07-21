# Demo 04: Biobank Screen-TRex (Algorithm 1) — In-Memory (Python)

## Purpose

Illustrates the **Biobank Screen-TRex workflow** ("Algorithm 1"), which screens
*many phenotypes* against a single shared design matrix X. Per phenotype the
algorithm tries **Screen-TRex Ordinary** first; if the estimated FDR falls
outside the acceptable window `[0.05, 0.15]` it retries with **Screen-TRex
Bootstrap-CI**; if that also lands outside the window it falls back to the
**classical T-Rex selector** at a target FDR of `0.10`. The interesting quantity
is the *routing behaviour*: which method handles which phenotype (the **Usage %**),
and how well each performs on its subset. This demo keeps the dummy matrices in
memory; Demo 05 runs the identical study memory-mapped. Python port of
[`cpp/trex_selector_methods/trex_screening/demo_trex_scr_04_mc_sim_biobank`](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_04_mc_sim_biobank/README.md).

## DGP / Data

`dgp_iid_snr` — i.i.d. Gaussian design, `n = 300`, `p = 1000`, `beta_j = 1`,
random support drawn per trial, SNR sweep `{0.1, 0.5, 1.0, 2.0, 5.0}`.

- **Part 1** — single phenotype, `s = 10`.
- **Part 2** — `q = 5` phenotypes sharing **one design matrix** X per trial
  (fixed `X_seed`), each with its own random support (`s = 5`) and noise draw —
  the biobank setting.

## Methods compared

The three routing outcomes of Algorithm 1 (all `ScreenTRexMethod.TREX`), keyed
in each result's `method_used` and shown under the display names:

- **Screen-TRex Ordinary: TLARS** — `"Screen-TRex (ordinary)"`.
- **Screen-TRex Bootstrap-CI: TLARS** — `"Screen-TRex (bootstrap-CI)"`.
- **T-Rex fallback: TLARS** — `"T-Rex (fallback)"`, at `target_FDR_trex = 0.10`.

FDR and TPR are **conditional on routing** (averaged only over the phenotypes a
method actually handled); the **Usage %** row is unconditional. A metric row for
a method with 0 % usage is an average over an empty set, not a performance
statement.

## What it computes

Each part sweeps the SNR grid; per SNR it runs `run_mc_biobank` over the trials
(`_biobank_single_worker` / `_biobank_multi_worker`) and aggregates with
`accumulate_snr`. Control: `K = 20`, window `[0.05, 0.15]`,
`target_FDR_trex = 0.10`. `num_MC = 60` (single) / `20` (multi) and
`R_boot = 500` are committed defaults — **downscaled** from the cpp demo's
`num_MC = 200` / `R_boot = 1000`, since the biobank workflow (bootstrap band +
classical T-Rex fallback) is expensive. Override with `SCR_NUM_MC` /
`SCR_NUM_MC_MULTI`. Trials run in a `ProcessPoolExecutor` (`SCR_NUM_WORKERS`,
default 6).

## Output

Written to `simulation_results/data/`:

- `scr_biobank_snr_n300_p1000_s10.txt` / `.csv` — Part 1.
- `scr_biobank_multi_n300_p1000_q5_s5.txt` / `.csv` — Part 2.

Each holds per-method Usage %, FDR, TPR, and estimated FDR (tidy CSV:
`method, metric, snr, value`, with a `Usage` fraction metric).

## Running

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_04_mc_sim_biobank/demo_trex_scr_04_mc_sim_biobank.py
```

## Interpretation

- **Part 1**: the routing does what it is designed to do — at very low SNR
  essentially every phenotype ends at the classical fallback; as SNR grows the
  traffic shifts to the screening routes, whose conditional error control is
  good once signal is present.
- **Part 2**: with only `s = 5` active variables per phenotype the fallback
  dominates, and the rarely-used screening routes are poorly calibrated (a
  single false selection moves the coarse FDP a long way) — the cautionary
  result of the pair.

**Last updated**: 2026-07-21
