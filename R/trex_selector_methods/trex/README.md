# Classical T-Rex Selector: R Demonstration Suite

## Purpose

R example programs for the classical **Terminating Random Experiments Selector
(T-Rex Selector)**, built on the **TRexSelectorNeo** R6 API:

```r
library(TRexSelectorNeo)

ctrl     <- trex_control(solver = "TLARS", K = 20L)
selector <- TRexSelector$new(X, y, tFDR = 0.1, control = ctrl)
selector$select()
selector$selected_indices          # 1-based selected support
compute_fdp(selector$selected_indices, true_support)
compute_tpp(selector$selected_indices, true_support)
```

The suite mirrors the C++ demos in
[cpp/trex_selector_methods/trex/](../../../cpp/trex_selector_methods/trex/README.md)
(see that README for the statistical model, notation, and FDR/TPR
definitions) and covers single runs, Monte Carlo studies with fixed and
variable support, l-loop strategy comparisons, memory-mapped workflows, and a
scalability benchmark.

---

## Demos

| # | File | Description | cpp counterpart |
|---|---|---|---|
| 00 | [demo_trex_00_gauss_data.R](demo_trex_00_gauss_data.R) | Runs the classical T-Rex (TLARS) on the package-shipped toy `Gauss_data` (n=50, p=100, s=3); reports calibration details, selected support, TPP and FDP | — (R-only) |
| 01 | [demo_trex_01_single_run.R](demo_trex_01_single_run.R) | Single run in a low- and a high-dimensional setting with a fixed true support (TLARS); illustrates reporting of the internal statistics `phi_prime`, `phi_mat`, `fdp_hat_mat`, `r_mat`, `voting_grid` | `demo_trex_01_single_run.cpp` |
| 02 | [demo_trex_02_mc_sim_fixed_support.R](demo_trex_02_mc_sim_fixed_support.R) | MC simulation with fixed support (drawn once with seed 24): 14 solvers x SNR sweep {0.1, 0.5, 1.0, 2.0, 5.0}; averaged FDR/TPR per solver x SNR | `demo_trex_02_mc_sim_fixed_support.cpp` |
| 03 | [demo_trex_03_mc_sim_variable_support.R](demo_trex_03_mc_sim_variable_support.R) | MC simulation with support (and optionally coefficients) redrawn each trial: 14 solvers x 21 SNR values {0.1, 0.2, …, 2.0, 5.0}; averaged FDR / TPR / Avg L / Avg T | `demo_trex_03_mc_sim_variable_support.cpp` |
| 04 | [demo_trex_04_mc_sim_lloop_strategies.R](demo_trex_04_mc_sim_lloop_strategies.R) | MC comparison of the six L-loop strategies (STANDARD, HCONCAT, PERMUTATION, PERMUTATION_DIRECT, DIRECT, SKIPL) with the TLARS base solver; n=300, p=1000, s=10; 21 SNR values; random and block support scenarios | `demo_trex_04_mc_sim_lloop_strategies.cpp` |
| 05 | [demo_trex_05_mmap.R](demo_trex_05_mmap.R) | Memory-mapped usage patterns. Part A: in-memory X + `use_memory_mapping = TRUE` (internal D-mmap + solver serialization). Part B: fully memory-mapped X via `convert_to_memory_mapped()` / `mmap_matrix()` | `demo_trex_05_mmap.cpp` |
| 06 | [demo_trex_06_mc_sim_mmap.R](demo_trex_06_mc_sim_mmap.R) | MC verification of the D-mmap + solver-serialization pipeline (TLARS only, fixed support seed 24, SNR {0.1, 0.5, 1.0, 2.0, 5.0}, `tloop_max_stagnant_steps = 7`, `opt_threshold = 0.75`). Part A: in-memory X; Part B: per-trial temporary mmap-backed X/y with exception-safe cleanup | `demo_trex_06_mc_sim_mmap.cpp` |
| 07 | [demo_trex_07_scalability.R](demo_trex_07_scalability.R) | MC scalability benchmark: in-memory vs chunked, memory-mapped out-of-core execution over an exponentially increasing (n, p) grid targeting specific memory footprints; measures execution time and peak memory, saves a CSV | — (R-only; the C++ demo 07 is a WIP placeholder) |

---

## Shared utilities: `trex_sim_utils.R`

[trex_sim_utils.R](trex_sim_utils.R) mirrors the C++ `trex_sim_utils.hpp` and
is sourced by every demo in this folder:

- **Section 1 — DGP (in-memory)**: `dgp_gauss_snr()` — i.i.d. N(0,1)
  predictors with SNR-calibrated noise.
- **Section 2 — DGP (out-of-core / mmap)**: `dgp_chunked_gauss_snr()` —
  chunked out-of-core DGP for massive datasets.
- **Section 3 — print helpers**: `.print_single_run()`.
- **Section 4 — MC infrastructure**:
  - `SOLVERS_DEFAULT` — 14 solver descriptors (TLARS, TLASSO, TENET,
    TSTAGEWISE, TSTEPWISE, TOMP, TGP, TACGP, TMP, TAFS with rho_afs 0.3 and
    1.0, TNCGMP variants v0/v1, TOOLS), mirroring
    `make_default_solvers_to_test()` in C++;
  - `make_mmap_trex_ctrl()` — `trex_control` builder for memory-mapped
    scenarios; `execute_with_temp_mmap()` — RAII-like lifecycle manager for
    temporary mmap files;
  - `.make_trex_ctrl()` — builds a `trex_control()` from a solver descriptor;
  - `.run_mc_trex()` — parallel MC runner over the R6 API
    (`parallel::mclapply` on Unix, `parLapply` on Windows);
  - `.print_solver_table()`, `.save_mc_csv()`, `.save_and_print_mc()` —
    console table plus combined `.txt` / `.csv` output.

Some demos additionally source the shared
[../support_generators.R](../support_generators.R) and
[../simulation_utils.R](../simulation_utils.R).

---

## Control parameters

All selector configuration goes through `trex_control()`. Parameters used
across these demos include:

- `solver` (e.g. `"TLARS"`), plus solver-specific `lambda2`, `rho_afs`,
  `ncgmp_variant`;
- `K` (number of random experiments) and `max_dummy_multiplier`;
- `lloop_strategy` (`"STANDARD"`, `"HCONCAT"`, `"PERMUTATION"`,
  `"PERMUTATION_DIRECT"`, `"DIRECT"`, `"SKIPL"`);
- `use_max_T_stop`, `tloop_stagnation_stop`, `tloop_max_stagnant_steps`,
  `opt_threshold`;
- `use_memory_mapping` (activates D-mmap + solver serialization);
- `use_openmp`, `seed`.

The target FDR is not part of the control object; it is passed as `tFDR` to
`TRexSelector$new()`.

---

## Seed conventions

- **Per-trial DGP seed**: `.run_mc_trex()` derives `trial_seed = base_seed +
  mc - 1` and hands it to the data-generating factory, so trials are
  reproducible and independent.
- **Dummy generation**: inside the MC runner, `trex_control(seed = -1L)` is
  used — hardware entropy per trial, giving independent random experiments.
- **Random supports**: `make_support_random()` seeds with
  `seed + 500000L` (mirrors the C++ `random_support()`).
- **Fixed supports**: demos 02 and 06 draw the support once with seed 24;
  demo 01 uses the fixed C++ support {27, 149, 43, 128, 42, 4} (0-based),
  converted to 1-based R indices {28, 150, 44, 129, 43, 5}.

---

## Running and outputs

```bash
Rscript R/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support.R
```

Demos run from any working directory. The MC demos (02, 03, 04, 06, 07) write
into [simulation_results/](simulation_results/) next to the demos:

- `.txt` files — human-readable summaries,
- `.csv` files — tidy long-format tables for plotting and post-processing
  (e.g. `demo_trex_02_mc_sim_fixed_support_results_n300_p1000_stagnation_window_7.csv`).

Demos 00, 01, and 05 print to the console only.

---

**Last updated**: 2026-07-06
