# Demo 05: Biobank Screen-TRex (Algorithm 1) — Memory-Mapped (Python)

## Purpose

Repeats the **Biobank Screen-TRex workflow** ("Algorithm 1") of
[Demo 04](../demo_trex_scr_04_mc_sim_biobank/README.md) with **memory-mapped
dummy matrices** (`use_memory_mapping = True`), the storage mode intended for
biobank-scale problems where holding the dummy matrices in RAM is infeasible.
The algorithm is unchanged (Ordinary → Bootstrap-CI → classical T-Rex fallback,
acceptable window `[0.05, 0.15]`, `target_FDR_trex = 0.10`); the point is that
moving the dummy matrices to disk changes the memory footprint and not the
statistics. Python port of
[`cpp/trex_selector_methods/trex_screening/demo_trex_scr_05_mc_sim_biobank_mmap`](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_05_mc_sim_biobank_mmap/README.md).

## DGP / Data

Identical to Demo 04: `dgp_iid_snr`, `n = 300`, `p = 1000`, SNR sweep
`{0.1, 0.5, 1.0, 2.0, 5.0}`. **Part 1** single phenotype (`s = 10`); **Part 2**
`q = 5` phenotypes sharing one X per trial (`s = 5`).

## Methods compared

The same three routing outcomes as Demo 04 (**Screen-TRex Ordinary: TLARS**,
**Screen-TRex Bootstrap-CI: TLARS**, **T-Rex fallback: TLARS**), with FDR/TPR
conditional on routing and Usage % unconditional.

## Memory mapping

The only difference from Demo 04 is `use_memory_mapping = True` in the nested
`trex_ctrl`. As in the cpp source, the design matrix **X still lives in RAM**;
only the dummy matrices move to disk, and the library manages those temporary
files in its own scratch directories — this script creates no memory-mapped
files of its own and has nothing to clean up.

## What it computes

The same structure as Demo 04 (`run_mc_biobank` + `accumulate_snr`). Committed
defaults `num_MC = 40` (single) / `15` (multi) and `R_boot = 500` are
**downscaled** from the cpp demo's `num_MC = 200` / `R_boot = 1000` (override
with `SCR_NUM_MC` / `SCR_NUM_MC_MULTI`). Trials run in a `ProcessPoolExecutor`
(`SCR_NUM_WORKERS`, default 6); running each memory-mapped selector inside a
fresh worker process also keeps the mmap lifetimes fully isolated.

## Output

Written to `simulation_results/data/`:

- `scr_biobank_mmap_snr_n300_p1000_s10.txt` / `.csv` — Part 1.
- `scr_biobank_mmap_multi_n300_p1000_q5_s5.txt` / `.csv` — Part 2.

## Running

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_05_mc_sim_biobank_mmap/demo_trex_scr_05_mc_sim_biobank_mmap.py
```

## Interpretation

- Routing and conditional error control behave **as in the in-memory Demo 04**;
  residual differences are Monte Carlo noise (`seed = -1` per trial), not a
  discrepancy in the memory-mapped path.
- **Memory mapping changes the footprint, not the statistics.**

**Last updated**: 2026-07-21
