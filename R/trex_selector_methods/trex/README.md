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
variable support, l-loop strategy and dummy distribution comparisons,
memory-mapped workflows, and a scalability benchmark.

The folder layout mirrors the C++ suite: one subfolder per demo (each with its
own `README.md` and `simulation_results/`), plus the suite-level
`trex_sim_utils.R` helper.

```txt
trex/
  ├── README.md
  ├── trex_sim_utils.R                 <- suite-level MC / DGP helpers
  ├── demo_trex_01_single_run/
  │   ├── demo_trex_01_single_run.R
  │   ├── README.md
  │   └── simulation_results/
  ├── demo_trex_02_mc_sim_fixed_support/
  │   └── ...
  └── demo_trex_08_mc_sim_scalability/
```

---

## Demos

| # | Folder | Description | cpp counterpart |
|---|---|---|---|
| 01 | [demo_trex_01_single_run](demo_trex_01_single_run/) | Single run in a low- and a high-dimensional setting with a fixed true support (TLARS); illustrates reporting of the internal statistics `phi_prime`, `phi_mat`, `fdp_hat_mat`, `r_mat`, `voting_grid` | `demo_trex_01_single_run` |
| 02 | [demo_trex_02_mc_sim_fixed_support](demo_trex_02_mc_sim_fixed_support/) | MC simulation with fixed support (drawn once with seed 24): 14 solvers x SNR sweep {0.1, 0.5, 1.0, 2.0, 5.0}; averaged FDR/TPR per solver x SNR | `demo_trex_02_mc_sim_fixed_support` |
| 03 | [demo_trex_03_mc_sim_variable_support](demo_trex_03_mc_sim_variable_support/) | MC simulation with support (and optionally coefficients) redrawn each trial: 14 solvers x 21 SNR values {0.1, 0.2, …, 2.0, 5.0}; averaged FDR / TPR / Avg L / Avg T | `demo_trex_03_mc_sim_variable_support` |
| 04 | [demo_trex_04_mc_sim_lloop_strategies](demo_trex_04_mc_sim_lloop_strategies/) | MC comparison of the six L-loop strategies (STANDARD, HCONCAT, PERMUTATION, PERMUTATION_ONDEMAND, ONDEMAND, SKIPL) with the TLARS base solver; n=300, p=1000, s=10; 21 SNR values; random and block support scenarios | `demo_trex_04_mc_sim_lloop_strategies` |
| 05 | [demo_trex_05_mc_sim_dummy_distributions](demo_trex_05_mc_sim_dummy_distributions/) | MC comparison of 12 dummy distribution variants (Normal, Uniform, Rademacher, StudentT df 3/5, Laplace, Gumbel, Triangle, Logistic, Mammen, sparse Rademacher s=0.1, UniformSphere dim=5) across 3 base solvers (TLARS/TOMP/TAFS rho_afs=1.0) with the STANDARD L-loop strategy held fixed; n=300, p=1000, s=10; 21 SNR values; one result file pair per solver | `demo_trex_05_mc_sim_dummy_distributions` |
| 06 | [demo_trex_06_mmap](demo_trex_06_mmap/) | Memory-mapped usage patterns. Part A: in-memory X + `use_memory_mapping = TRUE` (internal D-mmap + solver serialization). Part B: fully memory-mapped X via `convert_to_memory_mapped()` / `mmap_matrix()` | `demo_trex_06_mmap` |
| 07 | [demo_trex_07_mc_sim_mmap](demo_trex_07_mc_sim_mmap/) | MC verification of the D-mmap + solver-serialization pipeline (TLARS only, fixed support seed 24, SNR {0.1, 0.5, 1.0, 2.0, 5.0}, `tloop_max_stagnant_steps = 7`, `opt_threshold = 0.75`). Part A: in-memory X; Part B: per-trial temporary mmap-backed X/y with exception-safe cleanup | `demo_trex_07_mc_sim_mmap` |
| 08 | [demo_trex_08_mc_sim_scalability](demo_trex_08_mc_sim_scalability/) | MC scalability benchmark: in-memory vs chunked, memory-mapped out-of-core execution over an exponentially increasing (n, p) grid targeting specific memory footprints; measures execution time and peak memory, saves a CSV | `demo_trex_08_mc_sim_scalability` (R-only content; the C++ demo 08 is a WIP placeholder) |

Each demo folder has its own `README.md`.

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

Because each demo now lives one level deeper (in its own subfolder), the
suite-level helper is sourced as `../trex_sim_utils.R`; some demos additionally
source the cross-suite shared
[../../support_generators.R](../support_generators.R) and
[../../simulation_utils.R](../simulation_utils.R) (at the
`trex_selector_methods/` level).

---

## Control parameters

All selector configuration goes through `trex_control()`. Parameters used
across these demos include:

- `solver` (e.g. `"TLARS"`), plus solver-specific `lambda2`, `rho_afs`,
  `ncgmp_variant`;
- `K` (number of random experiments) and `max_dummy_multiplier`;
- `lloop_strategy` (`"STANDARD"`, `"HCONCAT"`, `"PERMUTATION"`,
  `"PERMUTATION_ONDEMAND"`, `"ONDEMAND"`, `"SKIPL"`);
- `dummy_distribution` (built with `trex_dummy_distribution()`; e.g.
  `trex_dummy_distribution("StudentT", df = 3)` — see demo 05);
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
- **Fixed supports**: demos 02 and 07 draw the support once with seed 24;
  demo 01 uses the fixed C++ support {27, 149, 43, 128, 42, 4} (0-based),
  converted to 1-based R indices {28, 150, 44, 129, 43, 5}.

---

## Running and outputs

```bash
Rscript R/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support/demo_trex_02_mc_sim_fixed_support.R
```

Demos run from any working directory. Unlike the `trex_da` / `trex_gvs` demos,
these scripts do **not** take a worker-cores argument — MC parallelism uses a
hardcoded `num_cores <- 6L` near the top of each MC demo (edit it to match your
machine). The MC demos write into each demo's own
`simulation_results/data/` subfolder (the `simulation_results/` umbrella holds
generated artifacts only — `data/` for tables, `plots/` for figures, matching
the C++ suite convention):

- `.txt` files — human-readable summaries,
- `.csv` files — tidy long-format tables for plotting and post-processing
  (e.g. `demo_trex_02_mc_sim_fixed_support/simulation_results/data/demo_trex_02_mc_sim_fixed_support_results_n300_p1000_stagnation_window_7.csv`).

Demos 02–08 persist result tables (`.csv`/`.txt`, or the mmap selection/
benchmark CSVs); only demo 01 (single run) is console-only.

---

**Last updated**: 2026-07-11
