# Demo 04: Biobank Screen-TRex (In-Memory)

## Purpose

Demonstrate the **Biobank Screen-TRex** workflow ("Algorithm 1"): screening **many phenotypes against one shared design matrix** `X`, with adaptive per-phenotype routing between Screen-TRex Ordinary, Screen-TRex Bootstrap-CI, and a classical T-Rex fallback. In-memory variant — see Demo 05 for the memory-mapped counterpart. Uses `TRexBiobankScreeningSelector`; the script sources the shared `../trex_scr_sim_utils.R`.

## Data-generating process

All parts use `dgp_iid_snr` / `.make_response` from the shared utils (i.i.d. Gaussian `X`, SNR-calibrated noise, `beta_j = 1`), with `n = 300`, `p = 1000`.

- **Part 1** (single phenotype) — true support `{5, 28, 43, 150, 399}` (1-based; C++ 0-based `{4,27,42,149,398}`), `s = 5`, `SNR = 1.0`.
- **Part 2** (multi-phenotype) — `q = 5` phenotypes sharing one fixed-seed `X`; per-phenotype support sizes `{5,6,5,7,5}` on their own random supports; per-phenotype `SNR in {1.0, 2.0, 5.0, 10.0, 20.0}`.
- **Part 3** (MC, single phenotype) — random support `s = 10` per trial, `SNR in {0.1, 0.5, 1.0, 2.0, 5.0}`.
- **Part 4** (MC, multi-phenotype) — `q = 5`, support `s = 5` per phenotype, same SNR grid, one shared `X` per MC run.

Selected indices are **1-based** in R (the C++ sources are 0-based).

## Methods compared

Adaptive **Algorithm 1** routing: for each phenotype the biobank selector routes adaptively against the acceptable-FDR window set by `trex_biobank_control(lower_bound_FDR = 0.05, upper_bound_FDR = 0.15)`:

1. **Screen-TRex (ordinary)** — tried first; accepted if its `estimated_FDR` lands inside `[0.05, 0.15]`.
2. **Screen-TRex (bootstrap-CI)** — tried if Ordinary's estimated FDR falls outside the window.
3. **T-Rex (fallback)** — classical T-Rex at the selector's `tFDR = 0.10` target, used if Bootstrap-CI also fails.

The fraction of phenotypes routed to each method is the reported **Usage %**. `select()` returns **one record for a vector `Y`** or **a plain list of records for a matrix `Y`** (one per phenotype, `res[[i]]` having the same shape as the single-phenotype return). Each record is a named list with `phenotype_index`, `method_used`, `estimated_FDR`, `estimated_FDR_screen_ordinary`, `estimated_FDR_screen_bootstrap`, `used_fallback_trex`, and `selected_indices` (1-based) plus its `_screen_ordinary` / `_screen_bootstrap` variants. This matches the Python binding (one `BiobankScreenTRexResult` for a 1-D `y`, a list for a 2-D `Y`).

## Parts

- **Part 1** — single phenotype: basic Algorithm 1 functionality, prints routing diagnostics for one dataset.
- **Part 2** — 5 phenotypes sharing one `X`: per-phenotype routing table plus a usage summary.
- **Part 3** — MC SNR sweep, single phenotype: per-method Usage % / FDR / TPR / Est. FDR.
- **Part 4** — MC SNR sweep, multi-phenotype: same per-method metrics.

The `run_part_*` flags at the top toggle individual parts.

## Running

```bash
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_04_biobank_inmem/demo_trex_scr_04_biobank_inmem.R
```

Runs from any working directory (no CLI worker-core argument). MC counts are reduced from the C++ reference's `num_MC = 500` to `NUM_MC_SINGLE = 40` (Part 3) and `NUM_MC_MULTI = 15` (Part 4, so `q x 15 = 75` phenotype screenings per SNR point) to keep runtime to roughly one minute. Results are statistically comparable to the C++ demo, **not** bit-identical — the two use different RNG streams.

The MC parts (3, 4) write a `.txt` summary and a tidy long-format `.csv` (with per-method Usage % rows) into the local `simulation_results/` folder; the single- and multi-phenotype table parts (1, 2) print to the console only.

## Interpretation

- The single- and multi-phenotype (non-MC) parts are best used to sanity-check the routing logic: confirm each phenotype's `method_used` is consistent with its `estimated_FDR` falling inside or outside the `[0.05, 0.15]` window.
- The committed MC results show the **fallback dominating at low SNR** and the **screen methods taking over as SNR rises**. In the single-phenotype sweep (`d04_...s10_inmem.csv`) the T-Rex fallback carries 100% of routing at `SNR = 0.1` and drops to 0% by `SNR = 2`, while Bootstrap-CI usage climbs to 0.75 (SNR 2) and 0.875 (SNR 5) — at weak signal Ordinary's estimated FDR drifts outside the window, so the routing escalates.
- Compare conditional FDP/TPP (averaged only over trials routed to a given method) against the **unconditional** estimated FDR: since each is conditioned on a different, non-random subset of trials, read the FDP/TPP columns in the context of that method's Usage %, not in isolation. Where a method's usage is 0 its FDP/TPP rows are just zeros.
- The multi-phenotype sweep (`d04_...q5_s5_inmem.csv`) shows the **T-Rex fallback carrying nearly all routing across the whole SNR grid** (Usage 0.81-1.0), unlike the single-phenotype case that shifts toward Bootstrap-CI at high SNR. Since all 5 phenotypes share the same `X`, this systematic difference points to cross-phenotype interaction (e.g. shared dummy/D-matrix reuse) rather than a genuine per-phenotype statistical difference — and is a useful signal to keep an eye on when scaling up phenotype counts.

**Last updated**: 2026-07-08
