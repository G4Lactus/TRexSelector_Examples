# Demo 05: Biobank Screen-TRex (Memory-Mapped)

## Purpose

The **Biobank Screen-TRex** workflow ("Algorithm 1") — adaptive per-phenotype routing between Screen-TRex Ordinary, Screen-TRex Bootstrap-CI, and a classical T-Rex fallback — run with the design matrix `X` and the internal dummy matrices backed by **memory-mapped files** (`use_memory_mapping = TRUE`), so a biobank-scale `X` never has to reside fully in RAM. Statistically identical to Demo 04: the DGP and the Ordinary → Bootstrap-CI → T-Rex-fallback routing (and its per-method Usage %) are unchanged; only the storage location of `X` and the D-matrices differs. See Demo 04 for the full description of the result accessors and `method_used` labels — they are identical here. Uses `TRexBiobankScreeningSelector`; the script sources the shared `../trex_scr_sim_utils.R`. Mirrors `cpp/trex_selector_methods/trex_screening/demo_trex_scr_05_biobank_screen_trex_mmap/`.

## Data-generating process

All parts use `dgp_iid_snr` / `.make_response` from the shared utils (i.i.d. Gaussian `X`, SNR-calibrated noise, `beta_j = 1`), with `n = 300`, `p = 1000`.

- **Part A** (single phenotype) — true support `{5, 28, 43, 150, 399}` (1-based; C++ 0-based `{4,27,42,149,398}`), `s = 5`, `SNR = 1.0`, data seed 123.
- **Part B** (single + multi phenotype) — the single case reuses Part A's dataset; the multi case has `q = 5` phenotypes sharing one fixed-seed `X`, per-phenotype support sizes `{5,6,5,7,5}` on their own random supports, per-phenotype `SNR in {1.0, 2.0, 5.0, 10.0, 20.0}`.
- **Part C** (MC, single phenotype) — random support `s = 10` per trial, `SNR in {0.1, 0.5, 1.0, 2.0, 5.0}`, `NUM_MC = 30`.

Selected indices are **1-based** in R (the C++ sources are 0-based).

## Methods compared

For each phenotype the biobank selector routes adaptively against the acceptable-FDR window set by `trex_biobank_control(lower_bound_FDR = 0.05, upper_bound_FDR = 0.15)`:

1. **Screen-TRex (ordinary)** — tried first; accepted if its `estimated_FDR` lands inside `[0.05, 0.15]`.
2. **Screen-TRex (bootstrap-CI)** — tried if Ordinary's estimated FDR falls outside the window.
3. **T-Rex (fallback)** — classical T-Rex at the selector's `tFDR = 0.10` target, used if Bootstrap-CI also fails.

The fraction of phenotypes routed to each method is the reported **Usage %**. Shared controls come via `make_scr_trex_ctrl(use_memory_mapping = TRUE)` (`solver = "TLARS"`, `K = 20`), `trex_screen_control(R_boot = ...)`, and `trex_biobank_control()`. The single fully-mmap `X` in Parts B and C is built with `TRexSelectorNeo::convert_to_memory_mapped()` into a temporary `.dat` file, wrapped in an RAII-like `with_temp_mmap()` helper that `unlink()`s the file on exit (see Demo 02).

## Parts

- **Part A** — single phenotype, **in-memory** `X` with `use_memory_mapping = TRUE` (D-matrices mmap only): basic Algorithm 1 functionality with routing diagnostics for one dataset.
- **Part B** — **fully memory-mapped** `X` via `with_temp_mmap()`: first a single phenotype on the mmap `X`, then 5 phenotypes sharing one mmap `X` with a per-phenotype routing table plus a usage summary.
- **Part C** — MC SNR sweep, single phenotype, mmap pipeline enabled (fresh temp mmap per trial): per-method Usage % / FDR / TPR / Est. FDR.

The `run_part_*` flags at the top toggle individual parts.

## Running

```bash
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_05_biobank_mmap/demo_trex_scr_05_biobank_mmap.R
```

Runs from any working directory (no CLI worker-core argument). The MC count is reduced from the C++ reference's `num_MC = 500` to `NUM_MC = 30` to keep runtime to roughly one minute; results are statistically comparable to the C++ demo, **not** bit-identical — the two use different RNG streams. Part C writes a `.txt` summary and a tidy long-format `.csv` (with per-method Usage % rows) as `simulation_results/d05_biobank_screen_trex_mc_snr_n300_p1000_s10_mmap.{txt,csv}`; Parts A and B print to the console only.

## Interpretation

- Parts A and B are best used to sanity-check the routing logic under memory-mapped storage: confirm each phenotype's `method_used` is consistent with its `estimated_FDR` falling inside or outside the `[0.05, 0.15]` window, and that both the D-mmap-only (Part A) and fully-mmap-`X` (Part B) paths run end-to-end. This demo is primarily a **feasibility / parity check** — that adaptive routing does not depend on assumptions that break when `X` and the D-matrices are disk-backed across many phenotype × trial combinations.
- The committed 30-run single-phenotype sweep shows the **fallback dominating at low SNR and the screen methods taking over as SNR rises**: T-Rex fallback carries 100% of routing at `SNR = 0.1`, drops to 3.3% at `SNR = 1`, and 0% by `SNR = 2`, while Bootstrap-CI usage climbs from 20% (SNR 0.5) to 73.3% (SNR 5) and Ordinary settles around 27-37% at high SNR — at weak signal Ordinary's estimated FDR drifts outside the window, so the routing escalates.
- Read conditional FDP/TPP (averaged only over trials routed to a given method) in the context of that method's **Usage %**, not in isolation, since each is conditioned on a different, non-random subset of trials. Where a method's usage is 0 its FDR/TPR rows are just zeros (e.g. the screen methods at `SNR = 0.1`). At `SNR >= 2` the screen methods reach TPR ~1.0 with FDR in the ~0.04-0.07 range, and their estimated FDR tracks the realized FDR reasonably (Bootstrap-CI ~0.065-0.069 estimated vs. ~0.044-0.058 realized); the fallback's estimated FDR is pinned at the `tFDR = 0.10` target.
- Read alongside Demo 04: since the DGP, routing logic, and control parameters are identical, these Usage % / FDP / TPP / estimated-FDR figures are expected to closely match Demo 04's at the same SNR points — the only intended difference is where `X` and the D-matrices live.

**Last updated**: 2026-07-08
