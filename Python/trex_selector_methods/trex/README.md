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

---

## Demos

| # | Script | Description | C++ counterpart |
|---|---|---|---|
| 01 | [demo_trex_01_single_run.py](demo_trex_01_single_run.py) | Single run of `ts.TRexSelector(...).select()` with the TLARS solver in one low-dimensional (n=5000, p=1000) and one high-dimensional (n=150, p=300) setting; fixed 6-element support; prints `phi_prime`, `phi_mat`, `fdp_hat_mat`, `r_mat`, and `voting_grid` via `print_single_run_result()` (console output only). | [demo_trex_01_single_run](../../../cpp/trex_selector_methods/trex/demo_trex_01_single_run/) |
| 02 | [demo_trex_02_mc_sim_fixed_support.py](demo_trex_02_mc_sim_fixed_support.py) | Monte Carlo study (`run_mc_trex()`, num_MC=100) over all 14 solver descriptors in `SOLVERS_DEFAULT` and SNR in {0.1, 0.5, 1, 2, 5}; n=300, p=1000, s=10; fixed support drawn once with seed 24 and shared across all trials; reports averaged FDR/TPR per solver x SNR. | [demo_trex_02_mc_sim_fixed_support](../../../cpp/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support/) |
| 03 | [demo_trex_03_mc_sim_variable_support.py](demo_trex_03_mc_sim_variable_support.py) | Monte Carlo study (num_MC=100) over the same 14 solvers and a finer SNR sweep {0.1, 0.2, ..., 2.0, 5.0} (21 values); support (and optionally coefficients) redrawn each trial; additionally tracks average L and T per solver x SNR. | [demo_trex_03_mc_sim_variable_support](../../../cpp/trex_selector_methods/trex/demo_trex_03_mc_sim_variable_support/) |
| 04 | [demo_trex_04_mc_sim_lloop_strategies.py](demo_trex_04_mc_sim_lloop_strategies.py) | L-loop strategy comparison with TLARS as the fixed base solver (num_MC=25, 21 SNR values, random or block support). The script documents six strategies (STANDARD, HCONCAT, PERMUTATION, PERMUTATION_DIRECT, DIRECT, SKIPL); as currently shipped only the `SKIPL_20p` and `SKIPL_50p` entries are enabled in `L_STRATEGIES` — uncomment the others to compare the full set. | [demo_trex_04_mc_sim_lloop_strategies](../../../cpp/trex_selector_methods/trex/demo_trex_04_mc_sim_lloop_strategies/) |
| 05 | [demo_trex_05_mmap.py](demo_trex_05_mmap.py) | Single-run memory-mapped workflows (seed=58). Part A: in-memory X with `use_memory_mapping=True` (D-mmap + solver LARS-path serialization). Part B: fully memory-mapped pipeline — X written to a temp binary file via `numpy_to_memmap()` and passed to `TRexSelector` as the zero-copy `to_numpy()` view. Each part runs a low-dim (n=5000, p=1000) and a high-dim (n=1000, p=5000) case and saves a per-variable selection CSV. | [demo_trex_05_mmap](../../../cpp/trex_selector_methods/trex/demo_trex_05_mmap/) |
| 06 | [demo_trex_06_mc_sim_mmap.py](demo_trex_06_mc_sim_mmap.py) | Monte Carlo study under memory mapping (TLARS only, num_MC=500, fixed support seed 24, SNR in {0.1, 0.5, 1, 2, 5}). Part C: in-memory X + `use_memory_mapping=True` via `run_mc_trex()`. Part D: fully mmapped X per trial via `run_mc_trex_full_mmap()` (each worker writes X to a unique temp file and deletes it afterwards). Both parts should produce statistically equivalent FDR/TPR — that equivalence is the verification goal. | [demo_trex_06_mc_sim_mmap](../../../cpp/trex_selector_methods/trex/demo_trex_06_mc_sim_mmap/) |

---

## Helper Modules

| Module | Contents |
|---|---|
| [dgp_gauss_snr.py](dgp_gauss_snr.py) | `dgp_gauss_snr(n, p, support, amplitude, coefs, snr, seed)` — data-generating process for the i.i.d. Gaussian scenario. Draws X with i.i.d. N(0,1) entries and calibrates the noise standard deviation to a target SNR (`noise_sigma = sqrt(Var(X @ beta) / snr)`, ddof=1, matching the C++ and R implementations). Returns a dict with `X`, `y`, `beta`, `true_support`, and the scenario scalars. |
| [support_generators.py](support_generators.py) | `make_support_capped_spread(s, p, max_gap)` — deterministic evenly-spaced support; `make_support_random(s, p, seed)` — uniformly random support without replacement (internally uses `seed + 500_000`, see seed conventions below). |
| [trex_helpers.py](trex_helpers.py) | Single-run helpers: `print_single_run_result()` (formatted result block including TP/FP counts, FDP/TPP, `phi_prime`, `phi_mat`, `fdp_hat_mat`, `r_mat`, `voting_grid`) and `save_selection_csv()` (per-variable CSV with columns `variable_index`, `phi_prime`, `selected`, `true_positive`; written via pandas). |
| [trex_sim_common.py](trex_sim_common.py) | Shared MC infrastructure: `SOLVERS_DEFAULT` (14 solver descriptors), `compute_fdp()` / `compute_tpp()`, `_make_trex_ctrl_from_dict()` (flat dict → `TRexControlParameter`), the `ProcessPoolExecutor` workers `_trial_worker()` / `_trial_worker_full_mmap()`, the parallel runners `run_mc_trex()` / `run_mc_trex_full_mmap()`, and the output helpers `print_solver_table()` / `save_mc_results()`. |

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

The single-run demos (01 and 05) instead pass a fixed integer seed to both
the DGP and `TRexSelector` for reproducibility — this is intentional for
single runs; the MC demos use `seed=-1`.

---

## Control Parameters

The demos configure the selector through plain dicts, converted inside each
worker process by `_make_trex_ctrl_from_dict()` into a
`ts.TRexControlParameter`:

- **Core fields**: `K` (dummy-matrix experiments per T-loop iteration),
  `max_dummy_multiplier` (max L = multiplier x p), `use_max_T_stop`,
  `opt_threshold`, `tloop_stagnation_stop`, `tloop_max_stagnant_steps`,
  `parallel_rnd_experiments`, `use_memory_mapping`.
- **Solver selection**: `ctrl.solver_type` is set from a
  `ts.SolverTypeForTRex` enum. The string names in `SOLVERS_DEFAULT` map to
  the 12 solver types: TLARS, TLASSO, TENET, TSTEPWISE, TSTAGEWISE, TOMP,
  TGP, TACGP, TMP, TNCGMP, TOOLS, TAFS (14 descriptors in total — TAFS is
  tested with `rho_afs` 0.3 and 1.0, TNCGMP with variants 0 and 1).
- **L-loop strategy**: `ctrl.lloop_strategy` is a `ts.LLoopStrategy` enum:
  STANDARD, HCONCAT, PERMUTATION, PERMUTATION_DIRECT, DIRECT, SKIPL.
- **Solver hyperparameters**: `ctrl.solver_params`
  (`SolverHyperparameters`) carries `lambda2` (TENET/ridge penalty),
  `rho_afs` (TAFS), and `ncgmp_variant` (TNCGMP).

Only flat dicts cross the process boundary — enum-bearing control objects are
built inside each worker, never pickled.

---

## Running the Demos

No build step is needed; run any demo directly:

```bash
python demo_trex_01_single_run.py
python demo_trex_02_mc_sim_fixed_support.py
```

Each demo inserts its own directory into `sys.path`, so the helper modules
resolve regardless of the current working directory. The MC demos parallelize
trials with `ProcessPoolExecutor` (`_NUM_WORKERS = 6` by default; adjust the
module-level constants `_NUM_WORKERS` / `_NUM_MC` to match your machine).

### Output files

Results are written to `simulation_results/` next to the demo script
(created on first use):

- **Demo 01**: console output only (no files).
- **Demos 02–04, 06** (via `save_mc_results()`): a plain-text results table
  (`<stem>.txt`) and a tidy long-format CSV (`<stem>.csv`, columns
  `solver`, `metric`, `snr`, `value` with metrics FDR/TPR and, where tracked,
  AvgL/AvgT). Stems encode the scenario, e.g.
  `demo_trex_02_mc_sim_fixed_support_results_n300_p1000_stagnation_window_7`,
  `d03_trex_mmap_demo_c_n300_p1000_sw7`.
- **Demo 05** (via `save_selection_csv()`): per-variable selection CSVs,
  e.g. `d02_trex_mmap_demo_a_n5000_p1000.csv` and
  `d02_trex_mmap_demo_b_n1000_p5000.csv`.

All CSVs are written with pandas and can be plotted with any tool of choice
(matplotlib, ggplot2, ...).

---

## References

- C++ suite: [cpp/trex_selector_methods/trex/](../../../cpp/trex_selector_methods/trex/README.md)
  (statistical background, notation, example results)
- R suite: `R/trex_selector_methods/trex/`
- Main package: `TRexSelector/Python/TRexSelectorNeo` (TRexSelectorNeo v0.2.0)

---

**Last updated**: 2026-07-06
