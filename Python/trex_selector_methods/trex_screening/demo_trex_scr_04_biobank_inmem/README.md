# Demo 04: Biobank Screen-TRex (In-Memory)

## Purpose

Demonstrate the **Biobank Screen-TRex** workflow ("Algorithm 1"): screening **many phenotypes against one shared design matrix** `X`, with adaptive per-phenotype routing between Screen-TRex Ordinary, Screen-TRex Bootstrap-CI, and a classical T-Rex fallback. In-memory variant — see Demo 05 for the memory-mapped counterpart. Uses `ts.TRexBiobankScreeningSelector`; the script bootstraps `sys.path` to import the shared `trex_scr_common.py` and `trex_scr_bb_common.py` from the parent suite directory.

## Data-generating process

All parts use `dgp_iid_snr` / `_make_response` from the shared utils (i.i.d. Gaussian `X`, SNR-calibrated noise, `beta_j = 1`), with `n = 300`, `p = 1000`.

- **Part 1** (single phenotype) — true support `{4, 27, 42, 149, 398}` (0-based), `s = 5`, `SNR = 1.0`, data seed 123.
- **Part 2** (multi-phenotype) — `q = 5` phenotypes sharing one fixed-seed `X`; per-phenotype support sizes `{5,6,5,7,5}` on their own random supports; per-phenotype `SNR in {1.0, 2.0, 5.0, 10.0, 20.0}`.
- **Part 3** (MC, single phenotype) — random support `s = 10` per trial, `SNR in {0.1, 0.5, 1.0, 2.0, 5.0}`, `NUM_MC_SINGLE = 40`.
- **Part 4** (MC, multi-phenotype) — `q = 5`, support `s = 5` per phenotype, same SNR grid, one shared `X` per MC run, `NUM_MC_MULTI = 15`.

Selected indices are **0-based** in Python (matching the C++ sources).

## Methods compared

Adaptive **Algorithm 1** routing: for each phenotype the biobank selector routes adaptively against the acceptable-FDR window set by `make_biobank_ctrl(lower_bound_FDR=0.05, upper_bound_FDR=0.15)`:

1. **Screen-TRex (ordinary)** — tried first; accepted if its `estimated_FDR` lands inside `[0.05, 0.15]`.
2. **Screen-TRex (bootstrap-CI)** — tried if Ordinary's estimated FDR falls outside the window.
3. **T-Rex (fallback)** — classical T-Rex at the selector's `target_FDR_trex = 0.10` target, used if Bootstrap-CI also fails.

The fraction of phenotypes routed to each method is the reported **Usage %**. `ts.TRexBiobankScreeningSelector(X, Y, bio_ctrl=..., seed=-1, verbose=False).select()` returns **one `BiobankScreenTRexResult`** for a 1-D response `y`, or a **`list[BiobankScreenTRexResult]`** for a 2-D response matrix `Y` (one per phenotype column). Each result object exposes `method_used` (the strings `"Screen-TRex (ordinary)"` / `"Screen-TRex (bootstrap-CI)"` / `"T-Rex (fallback)"`), `estimated_FDR`, `estimated_FDR_screen_ordinary`, `estimated_FDR_screen_bootstrap`, `selected_indices` (0-based), `used_fallback_trex`, and `phenotype_index`. This differs from the R sibling, which returns a list with a `$statistics` data frame plus a `$selected_indices` list — the Python API returns result objects directly, mirroring the C++ `BiobankScreenTRexResult` vector.

The control is built by `make_biobank_ctrl(use_mmap=, R_boot=, lower_bound_FDR=, upper_bound_FDR=, target_FDR_trex=)` in the shared `trex_scr_bb_common.py`. It produces a `BiobankScreenTRexControl` that nests `trex_screen_ctrl` (a `ScreenTRexControlParameter`), which in turn nests `trex_ctrl` (`solver = "TLARS"`, `K = 20`). This demo runs with `use_mmap=False`.

## Parts

- **Part 1** — single phenotype: basic Algorithm 1 functionality, prints routing diagnostics for one dataset.
- **Part 2** — 5 phenotypes sharing one `X`: per-phenotype routing table plus a usage summary.
- **Part 3** — MC SNR sweep, single phenotype: per-method Usage % / FDR / TPR / Est. FDR.
- **Part 4** — MC SNR sweep, multi-phenotype: same per-method metrics.

The `RUN_PART_*` flags at the top of the script toggle individual parts.

## Running

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_04_biobank_inmem/demo_trex_scr_04_biobank_inmem.py
```

Runs from any working directory (the script resolves the shared modules via a `sys.path` bootstrap of the demo directory and its parent suite directory). MC counts are reduced from the C++ reference's `num_MC = 500` to `NUM_MC_SINGLE = 40` (Part 3) and `NUM_MC_MULTI = 15` (Part 4, so `q x 15 = 75` phenotype screenings per SNR point) to keep runtime short. Results are statistically comparable to the C++ demo, **not** bit-identical — the two use different RNG streams.

The MC parts (3, 4) write a `.txt` summary and a tidy long-format `.csv` (with per-method Usage % rows) into the local `simulation_results/` folder as `d04_biobank_screen_trex_mc_snr_n300_p1000_s10_inmem.{txt,csv}` and `d04_biobank_screen_trex_mc_multi_n300_p1000_q5_s5_inmem.{txt,csv}`; the single- and multi-phenotype table parts (1, 2) print to the console only. `simulation_results/` ships empty (just a `.gitkeep`) until the demo is run.

## Interpretation

- The single- and multi-phenotype (non-MC) parts are best used to sanity-check the routing logic: confirm each phenotype's `method_used` is consistent with its `estimated_FDR` falling inside or outside the `[0.05, 0.15]` window.
- In the MC sweeps, expect the **fallback to dominate at low SNR** and the **screen methods to take over as SNR rises** — at weak signal Ordinary's estimated FDR drifts outside the window, so the routing escalates to Bootstrap-CI and then the T-Rex fallback.
- Compare conditional FDP/TPP (averaged only over trials routed to a given method) against the **unconditional** estimated FDR: since each is conditioned on a different, non-random subset of trials, read the FDP/TPP columns in the context of that method's Usage %, not in isolation. Where a method's usage is 0 its FDP/TPP rows are just zeros.
- Comparing the single- vs. multi-phenotype MC results is informative: since all 5 phenotypes in the multi-phenotype scenario share the same `X`, any systematic difference from the single-phenotype case points to cross-phenotype interaction (e.g. shared dummy/D-matrix reuse) rather than a genuine per-phenotype statistical difference — a useful signal to watch when scaling up phenotype counts.

**Last updated**: 2026-07-08
