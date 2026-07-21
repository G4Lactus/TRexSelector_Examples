# Classical T-Rex Selector: Python Demonstration Suite

## Purpose

This folder contains the Python example scripts for the classical
**Terminating Random Experiments Selector (T-Rex Selector)**, built on the
`trex_selector_neo` package. The demos mirror the C++ suite in
[cpp/trex_selector_methods/trex/](../../../cpp/trex_selector_methods/trex/)
and the R suite in `R/trex_selector_methods/trex/`:

1. how to generate synthetic Gaussian regression data with SNR-calibrated noise,
2. how to run the T-Rex Selector (`trex_selector_neo.TRexSelector`) in
   single-run and Monte Carlo settings,
3. how to interpret the resulting variable selections (FDP/TPP) and
   simulation summaries (averaged FDR/TPR, average L and T).

If you are new to this folder, start with **Demo 01** for a basic single run,
then continue with **Demo 02** for a Monte Carlo study across several
signal-to-noise ratios. See the
[C++ trex README](../../../cpp/trex_selector_methods/trex/README.md) for the
statistical model, notation, and FDR/TPR definitions — they apply unchanged
here.

> **Indexing note**: All support sets and selected indices in these demos are
> **0-based** (matching C++). The R package is 1-based; e.g. the Demo 01
> support `{27, 149, 43, 128, 42, 4}` corresponds to R's
> `{28, 150, 44, 129, 43, 5}`.

The folder layout mirrors the C++ suite: one subfolder per demo (each with its
own `README.md` and `simulation_results/`), plus the suite-level helper modules.

```txt
trex/
  ├── README.md
  ├── dgp_gauss_snr.py        ┐
  ├── support_generators.py   │  suite-level shared modules
  ├── trex_helpers.py         │  (imported by every demo)
  ├── trex_sim_common.py      ┘
  ├── demo_trex_01_single_run/
  │   ├── demo_trex_01_single_run.py
  │   ├── README.md
  │   └── simulation_results/
  ├── demo_trex_02_mc_sim_fixed_support/
  │   └── ...
  └── demo_trex_08_mc_sim_scalability/
```

---

## Demos

| # | Folder | Description | C++ counterpart |
|---|---|---|---|
| 01 | [demo_trex_01_single_run](demo_trex_01_single_run/) | Single run of `tsn.TRexSelector(...).select()` with the TLARS solver in one high-dimensional (n=150, p=300) and one low-dimensional (n=5000, p=1000) setting (run in that order, as in C++); fixed 6-element support; prints selected indices, FDP/TPP, sparse `phi_prime`/`phi_mat`, `fdp_hat_mat`, `r_mat`, and `voting_grid` via `print_single_run_result()`, dual console + `d01_trex_basic_n{n}_p{p}.txt` output. | [demo_trex_01_single_run](../../../cpp/trex_selector_methods/trex/demo_trex_01_single_run/) |
| 02 | [demo_trex_02_mc_sim_fixed_support](demo_trex_02_mc_sim_fixed_support/) | Monte Carlo study (`run_mc_trex()`, num_MC=100 — Python downscale of the C++ 200) over all 14 solver descriptors in `SOLVERS_DEFAULT` and SNR in {0.1, 0.5, 1, 2, 5}; n=300, p=1000, s=10; fixed support drawn once with seed 24 and shared across all trials; reports averaged FDR/TPR per solver x SNR. | [demo_trex_02_mc_sim_fixed_support](../../../cpp/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support/) |
| 03 | [demo_trex_03_mc_sim_variable_support](demo_trex_03_mc_sim_variable_support/) | Monte Carlo study (num_MC=100 — Python downscale of the C++ 200) over the same 14 solvers and a finer SNR sweep {0.1, 0.2, ..., 2.0, 5.0} (21 values); support (and optionally coefficients) redrawn each trial; additionally tracks average L and T per solver x SNR (stagnation window 5). | [demo_trex_03_mc_sim_variable_support](../../../cpp/trex_selector_methods/trex/demo_trex_03_mc_sim_variable_support/) |
| 04 | [demo_trex_04_mc_sim_lloop_strategies](demo_trex_04_mc_sim_lloop_strategies/) | L-loop strategy comparison across 3 base solvers (TLARS/TOMP/TAFS rho_afs=0.3; one result file pair per solver) over all eight strategy rows — STANDARD, HCONCAT, PERMUTATION, PERMUTATION_ONDEMAND, ONDEMAND, and SKIPL at 5p/10p/20p (num_MC=25 — Python downscale of the C++ 200; 21 SNR values, random or block support). Tests whether the fixed-budget SKIPL FDR overshoot is LARS-specific. | [demo_trex_04_mc_sim_lloop_strategies](../../../cpp/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies/) |
| 05 | [demo_trex_05_mc_sim_dummy_distributions](demo_trex_05_mc_sim_dummy_distributions/) | Dummy distribution comparison across 3 base solvers (TLARS/TOMP/TAFS rho_afs=0.3) with the STANDARD L-loop strategy held fixed (num_MC=10 — Python downscale of the C++ 200; SNR grid {0.1, 0.2, 0.5, 0.6, 1, 2, 5} mirroring the current C++ grid; one result file pair per solver). Sweeps 12 distribution rows (Normal, Uniform, Rademacher, StudentT df 3/5, Laplace, Gumbel, Triangle, Logistic, Mammen, sparse Rademacher s=0.1, UniformSphere dim=5) passed as picklable spec dicts and resolved to `tsn.DummyDistribution` inside each worker. | [demo_trex_05_mc_sim_dummy_distributions](../../../cpp/trex_selector_methods/trex/demo_trex_05_mc_sim_dummy_distributions/) |
| 06 | [demo_trex_06_mmap](demo_trex_06_mmap/) | Single-run memory-mapped workflows (seed=58). Part A: in-memory X with `use_memory_mapping=True` (D-mmap + solver LARS-path serialization). Part B: fully memory-mapped pipeline — X written to a temp binary file via `numpy_to_memmap()` and passed to `tsn.TRexSelector` as the zero-copy `to_numpy()` view. Each part runs a low-dim (n=5000, p=1000) and a high-dim (n=1000, p=5000) case and writes a dual-output `.txt` (C++-style selection block) plus a per-variable selection CSV. | [demo_trex_06_mmap](../../../cpp/trex_selector_methods/trex/demo_trex_06_mmap/) |
| 07 | [demo_trex_07_mc_sim_mmap](demo_trex_07_mc_sim_mmap/) | Monte Carlo study under memory mapping (TLARS only, num_MC=200 as in C++, fixed support seed 24, SNR in {0.1, 0.5, 1, 2, 5}). Part C: in-memory X + `use_memory_mapping=True` via `run_mc_trex()`. Part D: fully mmapped X per trial via `run_mc_trex_full_mmap()` (each worker writes X to a unique temp file and deletes it afterwards). Both parts should produce statistically equivalent FDR/TPR — that equivalence is the verification goal. | [demo_trex_07_mc_sim_mmap](../../../cpp/trex_selector_methods/trex/demo_trex_07_mc_sim_mmap/) |
| 08 | [demo_trex_08_mc_sim_scalability](demo_trex_08_mc_sim_scalability/) | Scalability benchmark on the fully memory-mapped pipeline: scenarios A (300 x 1k), B (1k x 10k), C (5k x 100k), 2D sweep over 3 base solvers (TLARS/TOMP/TAFS) x SNR {0.5, 1, 2}, num_MC=10 sequential trials per grid point. X generated chunk-wise into an mmap backing file; twelve timing/memory/quality metric rows per solver; incremental save per completed scenario to `demo_trex_08_scalability_results.{txt,csv}`. | [demo_trex_08_mc_sim_scalability](../../../cpp/trex_selector_methods/trex/demo_trex_08_mc_sim_scalability/) |

---

## Helper Modules

| Module | Contents |
|---|---|
| [dgp_gauss_snr.py](dgp_gauss_snr.py) | `dgp_gauss_snr(n, p, support, amplitude, coefs, snr, seed)` — data-generating process for the i.i.d. Gaussian scenario. Draws X with i.i.d. N(0,1) entries and calibrates the noise standard deviation to a target SNR (`noise_sigma = sqrt(Var(X @ beta) / snr)`, ddof=1, matching the C++ and R implementations). Returns a dict with `X`, `y`, `beta`, `true_support`, and the scenario scalars. |
| [support_generators.py](support_generators.py) | `make_support_capped_spread(s, p, max_gap)` — deterministic evenly-spaced support; `make_support_random(s, p, seed)` — uniformly random support without replacement (internally uses `seed + 500_000`, see seed conventions below). |
| [trex_helpers.py](trex_helpers.py) | Single-run helpers mirroring `trex_sim_utils.hpp`: `make_dual_out()` (console + file PrintFn), the sparse printers `print_sparse_vector()` / `print_sparse_matrix()` (nonzero entries with indices, dense fallback), `print_selection_results()` (C++-style block: selected indices, FDP/TPP, sparse `phi_prime`/`phi_mat`, `fdp_hat_mat`, `r_mat`, `voting_grid`), `print_single_run_result()` (data/calibration header + the selection block), and `save_selection_csv()` (per-variable CSV with columns `variable_index`, `phi_prime`, `selected`, `true_positive`; written via pandas). |
| [trex_sim_common.py](trex_sim_common.py) | Shared MC infrastructure: `SOLVERS_DEFAULT` (14 solver descriptors), `compute_fdp()` / `compute_tpp()`, `_make_trex_ctrl_from_dict()` (flat dict → `TRexControlParameter`), `_make_dummy_distribution()` (spec dict → `tsn.DummyDistribution`, see demo 05), the `ProcessPoolExecutor` workers `_trial_worker()` / `_trial_worker_full_mmap()`, the parallel runners `run_mc_trex()` / `run_mc_trex_full_mmap()`, and the output helpers `print_solver_table()` / `save_mc_results()`. |

---

## Seed-Space Conventions

The MC demos keep four RNG streams disjoint (see the header of
[trex_sim_common.py](trex_sim_common.py)):

| Stream | Seed | Mechanism |
|---|---|---|
| DGP (X and noise) | `trial_seed` | `dgp_gauss_snr(seed=trial_seed)` |
| Support (variable-support demos) | `trial_seed + 500_000` | `make_support_random(seed=trial_seed)` adds the offset internally |
| Coefficients (`rnd_coef=True`) | `trial_seed + 600_000` | `np.random.default_rng(trial_seed + 600_000)` |
| TRexSelector dummies | `seed=-1` (hardware entropy) | Required for valid MC FDR estimates — a fixed integer seed would produce structurally correlated dummies across trials |

Trial seeds are derived as `trial_seed = base_seed + mc - 1` for
`mc = 1, ..., num_MC`, with `base_seed = 24 + snr_idx * 1000` per SNR grid
point (mirrors the R convention). Fixed-support demos draw the support once
with `make_support_random(s, p, seed=24)`.

The single-run demos (01 and 06) instead pass fixed integer seeds to the DGP
and `TRexSelector` for reproducibility (demo 01 mirrors C++ with data seed 58
and selector seed 42; demo 06 uses 58 for both) — this is intentional for
single runs; the MC demos use `seed=-1`.

---

## Control Parameters

The demos configure the selector through plain dicts, converted inside each
worker process by `_make_trex_ctrl_from_dict()` into a
`tsn.TRexControlParameter`:

- **Core fields**: `K` (random experiments per T-loop iteration),
  `max_dummy_multiplier` (max L = multiplier x p), `use_max_T_stop`,
  `opt_threshold`, `tloop_stagnation_stop`, `tloop_max_stagnant_steps`,
  `parallel_rnd_experiments`, `use_memory_mapping`.
- **Solver selection**: `ctrl.solver_type` is set from a
  `tsn.SolverTypeForTRex` enum. The string names in `SOLVERS_DEFAULT` map to
  the 12 solver types: TLARS, TLASSO, TENET, TSTEPWISE, TSTAGEWISE, TOMP,
  TGP, TACGP, TMP, TNCGMP, TOOLS, TAFS (14 descriptors in total — TAFS is
  tested with `rho_afs` 0.3 and 1.0, TNCGMP with variants 0 and 1).
- **L-loop strategy**: `ctrl.lloop_strategy` is a `tsn.LLoopStrategy` enum:
  STANDARD, HCONCAT, PERMUTATION, PERMUTATION_ONDEMAND, ONDEMAND, SKIPL.
- **Dummy distribution**: `ctrl.dummy_distribution` is a
  `tsn.DummyDistribution`, built from a picklable spec dict by
  `_make_dummy_distribution()` (see demo 05).
- **Solver hyperparameters**: `ctrl.solver_params`
  (`SolverHyperparameters`) carries `lambda2` (TENET/ridge penalty),
  `rho_afs` (TAFS), and `ncgmp_variant` (TNCGMP).

Only flat dicts cross the process boundary — enum-bearing control objects are
built inside each worker, never pickled.

---

## Running the Demos

No build step is needed; run any demo directly (repo venv at `.venv/`):

```bash
python Python/trex_selector_methods/trex/demo_trex_01_single_run/demo_trex_01_single_run.py
python Python/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support/demo_trex_02_mc_sim_fixed_support.py
```

Each demo inserts both its own directory **and its parent suite directory**
into `sys.path`, so the suite-level helper modules resolve regardless of the
current working directory (and are re-importable by the `spawn`-launched worker
processes). The MC demos parallelize trials with `ProcessPoolExecutor`
(`_NUM_WORKERS = 6` by default; adjust the module-level constants
`_NUM_WORKERS` / `_NUM_MC` to match your machine).

### Output files

Results are written to `simulation_results/data/` next to the demo script
(created on first use; the `simulation_results/` umbrella holds generated
artifacts only — `data/` for tables, `plots/` for figures, matching the C++
suite convention):

- **Demo 01**: dual console + file output, one `.txt` per scenario
  (`d01_trex_basic_n{n}_p{p}.txt`, C++-style selection block with sparse
  Phi_prime/Phi printing).
- **Demos 02–05, 07** (via `save_mc_results()`): a plain-text results table
  (`<stem>.txt`) and a tidy long-format CSV (`<stem>.csv`, columns
  `solver`, `metric`, `snr`, `value` with metrics FDR/TPR and, where tracked,
  AvgL/AvgT). Stems match the C++ result files, e.g.
  `demo_trex_02_mc_sim_fixed_support_trex_results_n300_p1000_stagnation_window_7`,
  `demo_trex_04_lloop_strategies_results_n300_p1000_tlars_random_support`,
  `d03_trex_mmap_demo_c_n300_p1000_sw7`.
- **Demo 06**: per configuration a dual-output `.txt` (C++-style selection
  block) plus a per-variable selection CSV (via `save_selection_csv()`),
  e.g. `d02_trex_mmap_demo_a_n5000_p1000.{txt,csv}` and
  `d02_trex_mmap_demo_b_n1000_p5000.{txt,csv}`.
- **Demo 08**: an aligned scenario x SNR metrics table
  (`demo_trex_08_scalability_results.txt`) and a tidy CSV with columns
  `scenario`, `solver`, `n`, `p`, `num_mc`, `snr`, `metric`, `value`,
  re-saved after each completed scenario.

All CSVs are written with pandas and can be plotted with any tool of choice
(matplotlib, ggplot2, ...).

---

## References

- C++ suite: [cpp/trex_selector_methods/trex/](../../../cpp/trex_selector_methods/trex/README.md)
  (statistical background, notation, example results)
- R suite: `R/trex_selector_methods/trex/`
- Main package: `TRexSelector/Python/TRexSelectorNeo` (TRexSelectorNeo v0.2.0)

---

**Last updated**: 2026-07-21
