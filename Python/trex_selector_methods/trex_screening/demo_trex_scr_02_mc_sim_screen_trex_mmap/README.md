# Demo 02: Screen-TRex with Memory-Mapped Dummy Matrices (Python)

## Purpose

Shows that the **memory-mapped Screen-TRex workflow** reproduces the in-memory
baseline while keeping the dummy (D) matrices on disk instead of in RAM. The
study is the same experiment as [Demo 01](../demo_trex_scr_01_mc_sim_screen_trex/README.md)
— same i.i.d. Gaussian design, same SNR sweep, same two thresholding rules — with
the single change `TRexControlParameter.use_memory_mapping = True`. The point is
*equivalence at a lower memory footprint*, not a different statistical regime.
Python port of
[`cpp/trex_selector_methods/trex_screening/demo_trex_scr_02_mc_sim_screen_trex_mmap`](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_02_mc_sim_screen_trex_mmap/README.md).

## DGP / Data

Identical to Demo 01: `dgp_iid_snr`, `n = 300`, `p = 1000`, `s = 10`, random
support drawn fresh per trial (`beta_j = 1`), SNR sweep
`{0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.

## Methods compared

The two `ScreenTRexMethod.TREX` rules — **Screen-TRex Ordinary: TLARS** and
**Screen-TRex Bootstrap-CI: TLARS** (`R_boot = 1000`, `ci_grid_step = 0.001`) —
with `use_memory_mapping = True` in the shared solver control (`solver = "TLARS"`,
`K = 20`).

## Memory mapping

As in the cpp source, the design matrix **X still lives in RAM**; only the dummy
matrices move to disk. The `use_memory_mapping = True` flag activates the
library's internal D-mmap + solver-serialization pipeline, which manages its own
temporary scratch files — this script creates no memory-mapped files of its own
and has nothing to clean up.

## What it computes

`demo_screen_trex_mmap_monte_carlo()` sweeps the two methods x the SNR grid,
`base_seed = 24 + snr_idx * 1000`, recording mean FDR, TPR, and estimated FDR.
`num_MC = 100` (committed default; cpp uses 200 — downscaled, override with
`SCR_NUM_MC`). The MC loop is parallelized with `ProcessPoolExecutor`
(`SCR_NUM_WORKERS`, default 6).

## Output

Written to `simulation_results/data/`:

- `scr_screen_trex_mmap_snr_n300_p1000.txt` / `.csv` — FDR, TPR, and estimated
  FDR per method and SNR level.

## Running

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_02_mc_sim_screen_trex_mmap/demo_trex_scr_02_mc_sim_screen_trex_mmap.py
```

## Interpretation

- **Memory mapping changes the memory footprint, not the statistical
  behaviour**: the mmap curves track the in-memory Demo 01 curves closely.
- **Exact equality is not expected** — both demos use nondeterministic selector
  seeds (`seed = -1`), so residual gaps are Monte Carlo noise, not a discrepancy
  in the memory-mapped path.
- The baseline conclusions carry over unchanged (uncontrolled at low SNR,
  reliable from `SNR >= 0.5`, Bootstrap-CI the conservative rule).

**Last updated**: 2026-07-21
