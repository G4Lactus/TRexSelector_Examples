# Screen-TRex Selector: R Demonstration Suite

## Purpose

R example programs for **Screen-TRex**, the T-Rex Selector extended with a
variable-**screening** step for ultra-high-dimensional data. Screen-TRex is
built for the regime where `p` is so large that running the classical T-Rex
selector directly is impractical (genome- or biobank-scale designs): it first
thresholds a dummy-based voting statistic `Phi_j` to screen down to a small
candidate set, using either the **Ordinary** majority-vote rule
(`Phi_j > 0.5`) or a **Bootstrap-CI** rule, plus **dependency-aware (DA)**
variants that recover power under correlated designs.

```r
library(TRexSelectorNeo)

# Ordinary / DA screening
sel <- TRexScreeningSelector$new(
  X, y, seed = 42L, verbose = FALSE,
  screen_control = trex_screen_control(trex_method = "TREX"),  # or a DA variant
  control        = trex_control(solver = "TLARS", K = 20L)
)
sel$select()
sel$selected_indices        # 1-based screened support
sel$estimated_FDR           # procedure's self-reported FDR
sel$estimated_correlation   # estimated rho (DA variants)
compute_fdp(sel$selected_indices, true_support)
compute_tpp(sel$selected_indices, true_support)

# Biobank / Algorithm 1 (many phenotypes, one shared X)
bio <- TRexBiobankScreeningSelector$new(
  X, Y, tFDR = 0.10, seed = 42L,
  biobank_control = trex_biobank_control(lower_bound_FDR = 0.05, upper_bound_FDR = 0.15),
  screen_control  = trex_screen_control(),
  control         = trex_control()
)
res <- bio$select()     # vector Y -> one record; matrix Y -> list of records
res[[1]]$method_used         # routing decision for phenotype 1
res[[1]]$selected_indices    # phenotype 1's 1-based selected indices
res[[1]]$estimated_FDR
```

The suite mirrors the C++ demos in
[cpp/trex_selector_methods/trex_screening/](../../../cpp/trex_selector_methods/trex_screening/README.md)
(see that README for the statistical model, notation, `Phi`, DA screening, and
Algorithm 1). Support sets and selected indices are **1-based** in R
(the C++ sources are 0-based).

The `trex_method` strings accepted by `trex_screen_control()` are:
`"TREX"` (ordinary voting-based screening), `"TREX_DA_AR1"`, `"TREX_DA_EQUI"`,
and `"TREX_DA_BLOCK_EQUI"` (dependency-aware variants). An unknown string
raises `Unknown ScreenTRexMethod: ...` at `select()` time.

The folder layout mirrors the C++ suite one-to-one: one subfolder per demo,
plus the shared `trex_scr_sim_utils.R` helper at the suite level.

```txt
trex_screening/
  ├── README.md
  ├── trex_scr_sim_utils.R            <- shared DGPs, MC runner, table output
  ├── demo_trex_scr_01_screen/
  │   ├── demo_trex_scr_01_screen.R
  │   ├── README.md
  │   └── simulation_results/
  ├── demo_trex_scr_02_screen_mmap/
  │   └── ...
  ├── demo_trex_scr_03_screen_correlated/
  ├── demo_trex_scr_04_biobank_inmem/
  ├── demo_trex_scr_05_biobank_mmap/
  └── demo_trex_scr_06_screen_solvers/
```

---

## Demos

Each demo folder has its own `README.md`.

| # | Folder | Description | cpp counterpart | Key APIs |
|---|---|---|---|---|
| 01 | [demo_trex_scr_01_screen/](demo_trex_scr_01_screen/) | In-memory Ordinary vs. Bootstrap-CI screening: two single runs (fixed support, SNR 5) plus an 8-point SNR-sweep Monte Carlo reporting FDR / TPR / Estimated FDR | [demo_trex_scr_01_screen_trex/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_01_screen_trex/) | `TRexScreeningSelector`, `trex_screen_control(use_bootstrap_CI=)` |
| 02 | [demo_trex_scr_02_screen_mmap/](demo_trex_scr_02_screen_mmap/) | Same screening scenario with the memory-mapped pipeline: Part A in-memory X + `use_memory_mapping=TRUE` (D-mmap); Part B fully memory-mapped X via `convert_to_memory_mapped()`; Part C mmap MC SNR sweep | [demo_trex_scr_02_screen_trex_mmap/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_02_screen_trex_mmap/) | `convert_to_memory_mapped`, `trex_control(use_memory_mapping=)` |
| 03 | [demo_trex_scr_03_screen_correlated/](demo_trex_scr_03_screen_correlated/) | Ordinary vs. DA screening on correlated designs: 4 single runs (Ordinary/DA-AR1/DA-EQUI/DA-BLOCK-EQUI) + 5 Monte Carlo sweeps (AR(1)/equi SNR sweeps and AR(1)/equi/block-equi rho sweeps) | [demo_trex_scr_03_screen_trex_correlated/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_03_screen_trex_correlated/) | `trex_method = "TREX_DA_*"`, `estimated_correlation` |
| 04 | [demo_trex_scr_04_biobank_inmem/](demo_trex_scr_04_biobank_inmem/) | Biobank "Algorithm 1" (in-memory): single phenotype, 5-phenotype routing table + summary, and MC SNR sweeps (single and multi-phenotype) reporting per-method Usage % / FDR / TPR / Est. FDR | [demo_trex_scr_04_biobank_screen_trex_inmem/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_04_biobank_screen_trex_inmem/) | `TRexBiobankScreeningSelector`, `trex_biobank_control` |
| 05 | [demo_trex_scr_05_biobank_mmap/](demo_trex_scr_05_biobank_mmap/) | Same biobank workflow with the memory-mapped pipeline: Part A in-memory X + D-mmap; Part B single + multi phenotype on a fully memory-mapped X; Part C mmap MC SNR sweep | [demo_trex_scr_05_biobank_screen_trex_mmap/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_05_biobank_screen_trex_mmap/) | `TRexBiobankScreeningSelector` + mmap X |
| 06 | [demo_trex_scr_06_screen_solvers/](demo_trex_scr_06_screen_solvers/) | Solver-backend comparison under screening: 8-point SNR sweep for TLARS / TAFS(rho_afs=0.3) / TOMP, Part 1 Ordinary only (3 series), Part 2 x {Ordinary, Bootstrap-CI} (6 series) | [demo_trex_scr_06_screen_trex_solvers/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_trex_solvers/) | `trex_control(solver=, rho_afs=)` |

---

## Shared utilities: `trex_scr_sim_utils.R`

[trex_scr_sim_utils.R](trex_scr_sim_utils.R) mirrors the two C++ demo-internal
headers (`demo_trex_scr_common.hpp`, `demo_trex_scr_bb_common.hpp`) and is
sourced by every demo:

- **Section 1 — DGPs**: `dgp_iid_snr()` (i.i.d. Gaussian), `dgp_ar1()`,
  `dgp_equi()`, `dgp_block_equi()` (AR(1) / equi / block-equi correlated
  designs), and `make_beta()` (evenly-spaced sparse coefficients). All use
  SNR-calibrated noise, `sigma = sqrt(||X beta||^2 / (n * SNR))`.
- **Section 2 — control factories**: `make_scr_trex_ctrl()`,
  `default_scr_methods()`.
- **Section 3 — single-run printer**: `.print_scr_result()`.
- **Section 4 — Monte Carlo**: `run_mc_screen()` (per-grid-point MC runner) and
  `save_and_print_scr_mc()` (console table + `.txt` + tidy `.csv`).
- **Section 5 — Biobank**: `.print_biobank_single()`, `.print_biobank_table()`,
  `.print_biobank_summary()`, and `save_and_print_biobank_mc()` (with a
  per-method Usage (%) row).

---

## Control parameters

- `trex_screen_control(use_bootstrap_CI, R_boot, ci_grid_step, trex_method)` —
  the screening rule and variant. **Note**: unlike the C++
  `make_screen_control()`, the R binding does **not** expose `rho_thr_DA` or
  `n_blocks`; the DA variants (including `TREX_DA_BLOCK_EQUI`) use the library's
  internal defaults for those, so the block DGP's `n_blocks` is not forwarded to
  the screener.
- `trex_control(solver, K, rho_afs, use_memory_mapping, ...)` — the underlying
  T-Rex solver configuration (same object as the classical `trex/` demos).
- `trex_biobank_control(lower_bound_FDR, upper_bound_FDR)` — the acceptable-FDR
  window for Algorithm 1's Ordinary -> Bootstrap-CI -> T-Rex-fallback routing.
  The T-Rex fallback target is the selector's `tFDR` argument.

The biobank `select()` returns **one record for a vector `Y`** (single
phenotype) or **a plain list of records for a matrix `Y`** (one per phenotype,
`res[[i]]` having the same shape as the single-phenotype return). Each record is
a named list with fields `phenotype_index`, `estimated_FDR`, `method_used`
(`"Screen-TRex (ordinary)"`, `"Screen-TRex (bootstrap-CI)"`,
`"T-Rex (fallback)"`), `used_fallback_trex`, `estimated_FDR_screen_ordinary`,
`estimated_FDR_screen_bootstrap`, `selected_indices` (1-based), and
`selected_indices_screen_ordinary` / `selected_indices_screen_bootstrap`. This
matches the Python binding (one `BiobankScreenTRexResult` for a 1-D `y`, or a
list of them for a 2-D `Y`). A `$statistics`-style summary table is a one-line
projection:
`do.call(rbind, lapply(res, function(r) as.data.frame(r[c("phenotype_index", "estimated_FDR", "method_used", "used_fallback_trex")])))`.

---

## Monte Carlo counts

The C++ demos use `num_MC = 500`. To keep each R demo to roughly one to two
minutes single-threaded, the MC counts here are reduced (60 for demo 01, 40 for
demo 02, 10 for demo 03's five sweeps, 40/15 for demo 04, 30 for demo 05, 20 for
demo 06) — each demo header states its own count. Results are statistically
comparable to the C++ reference, **not** bit-identical: the two use different
RNG streams, so seeds are mirrored as labels only.

---

## Running and outputs

```bash
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_01_screen/demo_trex_scr_01_screen.R
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_02_screen_mmap/demo_trex_scr_02_screen_mmap.R
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_03_screen_correlated/demo_trex_scr_03_screen_correlated.R
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_04_biobank_inmem/demo_trex_scr_04_biobank_inmem.R
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_05_biobank_mmap/demo_trex_scr_05_biobank_mmap.R
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_solvers/demo_trex_scr_06_screen_solvers.R
```

Demos run from any working directory (each resolves the shared
`../trex_scr_sim_utils.R` relative to its own location). The Monte Carlo parts
write into each demo's own `simulation_results/` subfolder:

- `.txt` — human-readable summary tables,
- `.csv` — tidy long-format tables (`method, metric, SNR|rho, value`) for
  plotting and post-processing.

The single-run parts (and the biobank single/multi-phenotype tables) print to
the console only.

---

**Last updated**: 2026-07-08
