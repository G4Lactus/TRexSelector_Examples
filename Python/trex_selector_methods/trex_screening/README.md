# Screen-TRex Selector: Python Demonstration Suite

## Purpose

Python example programs for **Screen-TRex**, the T-Rex Selector extended with a
variable-**screening** step for ultra-high-dimensional data. Screen-TRex is
built for the regime where `p` is so large that running the classical T-Rex
selector directly is impractical (genome- or biobank-scale designs): it first
thresholds a dummy-based voting statistic `Phi_j` to screen down to a small
candidate set, using either the **Ordinary** majority-vote rule
(`Phi_j > 0.5`) or a **Bootstrap-CI** rule, plus **dependency-aware (DA)**
variants that recover power under correlated designs. Unlike the classical
T-Rex selector, screening has **no target-FDR parameter**; it reports its own
*estimated FDR* alongside the selection, and whether that self-estimate can be
trusted depends sharply on the correlation structure (see Demos 01 and 03).

All six demos are **Monte Carlo studies** that parallelize their trials with
`concurrent.futures.ProcessPoolExecutor` and write per-method sweep tables
(`.txt` + tidy `.csv`) into each demo's own `simulation_results/data/` folder.

```python
import trex_selector_neo as tsn

# Ordinary / DA screening
sc = tsn.ScreenTRexControlParameter()
sc.trex_method = tsn.ScreenTRexMethod.TREX          # or a TREX_DA_* variant
tc = tsn.TRexControlParameter(); tc.K = 20
sel = tsn.TRexScreeningSelector(X, y, screen_control=sc, trex_control=tc,
                                seed=-1, verbose=False)
res = sel.select()
sel.selected_indices          # 0-based screened support (numpy array)
res.estimated_FDR             # procedure's self-reported FDR
res.estimated_correlation     # estimated rho (DA variants)

from trex_selector_neo.utils import compute_fdp, compute_tpp
compute_fdp(list(sel.selected_indices), true_support)
compute_tpp(list(sel.selected_indices), true_support)

# Biobank / Algorithm 1 (many phenotypes, one shared X)
bc = tsn.BiobankScreenTRexControl()
bc.lower_bound_FDR = 0.05; bc.upper_bound_FDR = 0.15; bc.target_FDR_trex = 0.10
bio = tsn.TRexBiobankScreeningSelector(X, Y, bio_ctrl=bc, seed=-1, verbose=False)
results = bio.select()        # BiobankScreenTRexResult, or a list for m-D Y
results[0].method_used        # routing decision per phenotype
results[0].selected_indices   # per-phenotype 0-based selected indices
```

The suite mirrors the C++ demos in
[cpp/trex_selector_methods/trex_screening/](../../../cpp/trex_selector_methods/trex_screening/README.md)
one-to-one (see that README for the statistical model, notation, `Phi`, DA
screening, Algorithm 1, and the annotated result discussion). Support sets and
selected indices are **0-based** in Python (matching the C++ sources; the R
suite is 1-based).

The `ScreenTRexMethod` values accepted by `ScreenTRexControlParameter.trex_method`
are `tsn.ScreenTRexMethod.TREX` (ordinary voting-based screening), `.TREX_DA_AR1`,
`.TREX_DA_EQUI`, and `.TREX_DA_BLOCK_EQUI` (dependency-aware variants).

The folder layout mirrors the C++ suite one-to-one: one subfolder per demo,
plus the two shared helper modules at the suite level.

```txt
trex_screening/
  ├── README.md
  ├── trex_scr_common.py              <- DGPs, control builders, MC runner, table output
  ├── trex_scr_bb_common.py           <- Biobank (Algorithm 1) MC runner + table output
  ├── demo_trex_scr_01_mc_sim_screen_trex/
  │   ├── demo_trex_scr_01_mc_sim_screen_trex.py
  │   ├── README.md
  │   └── simulation_results/
  ├── demo_trex_scr_02_mc_sim_screen_trex_mmap/
  │   └── ...
  ├── demo_trex_scr_03_mc_sim_correlated/
  ├── demo_trex_scr_04_mc_sim_biobank/
  ├── demo_trex_scr_05_mc_sim_biobank_mmap/
  └── demo_trex_scr_06_mc_sim_solvers/
```

---

## Demos

Each demo folder has its own `README.md`.

| # | Folder | Description | cpp counterpart | Key APIs |
|---|---|---|---|---|
| 01 | [demo_trex_scr_01_mc_sim_screen_trex/](demo_trex_scr_01_mc_sim_screen_trex/) | In-memory Ordinary vs. Bootstrap-CI screening: an 8-point SNR-sweep Monte Carlo reporting FDR / TPR / Estimated FDR | [demo_trex_scr_01_mc_sim_screen_trex/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_01_mc_sim_screen_trex/) | `TRexScreeningSelector`, `ScreenTRexControlParameter.use_bootstrap_CI` |
| 02 | [demo_trex_scr_02_mc_sim_screen_trex_mmap/](demo_trex_scr_02_mc_sim_screen_trex_mmap/) | The same SNR sweep with `use_memory_mapping=True` (dummy matrices on disk); equivalence to Demo 01 at a lower memory footprint | [demo_trex_scr_02_mc_sim_screen_trex_mmap/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_02_mc_sim_screen_trex_mmap/) | `TRexControlParameter.use_memory_mapping` |
| 03 | [demo_trex_scr_03_mc_sim_correlated/](demo_trex_scr_03_mc_sim_correlated/) | Ordinary/Bootstrap vs. DA screening on correlated designs: 5 MC sweeps (AR(1)/equi SNR sweeps and AR(1)/equi/block-equi rho sweeps) | [demo_trex_scr_03_mc_sim_correlated/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_03_mc_sim_correlated/) | `ScreenTRexMethod.TREX_DA_*`, `estimated_correlation` |
| 04 | [demo_trex_scr_04_mc_sim_biobank/](demo_trex_scr_04_mc_sim_biobank/) | Biobank "Algorithm 1" (in-memory): MC SNR sweeps (single and multi-phenotype) reporting per-method Usage % / FDR / TPR / Est. FDR | [demo_trex_scr_04_mc_sim_biobank/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_04_mc_sim_biobank/) | `TRexBiobankScreeningSelector`, `BiobankScreenTRexControl` |
| 05 | [demo_trex_scr_05_mc_sim_biobank_mmap/](demo_trex_scr_05_mc_sim_biobank_mmap/) | The same biobank MC sweeps with `use_memory_mapping=True` | [demo_trex_scr_05_mc_sim_biobank_mmap/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_05_mc_sim_biobank_mmap/) | `TRexBiobankScreeningSelector` + D-mmap |
| 06 | [demo_trex_scr_06_mc_sim_solvers/](demo_trex_scr_06_mc_sim_solvers/) | Solver-backend comparison under screening: SNR sweep for TLARS / TAFS(rho_afs=0.3) / TOMP, Part 1 Ordinary only (3 series), Part 2 x {Ordinary, Bootstrap-CI} (6 series) | [demo_trex_scr_06_mc_sim_solvers/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_mc_sim_solvers/) | `TRexControlParameter.solver_type`, `solver_params.rho_afs` |

---

## Shared modules

Both live at the suite level and mirror the C++ demo-internal headers
(`trex_screening_dgps.hpp`, `trex_screening_simulation_utils.hpp`). Each demo
resolves them via a `sys.path` bootstrap that inserts both its own directory and
the parent (suite) directory, so the shared modules import cleanly from the
nested demo location and — because macOS spawns worker processes — from every
`ProcessPoolExecutor` child.

- [trex_scr_common.py](trex_scr_common.py) — DGPs (`dgp_iid_snr`, `dgp_ar1`,
  `dgp_equi`, `dgp_block_equi`, `make_beta`), flat-dict control builders
  (`_build_trex_ctrl`, `_build_screen_ctrl`, `default_scr_methods`), the
  module-level ProcessPoolExecutor worker `_screen_trial_worker`, the parallel
  Monte Carlo runner `run_mc_screen`, and the table writer
  `save_and_print_scr_mc` (console + `.txt` + tidy pandas `.csv`). All DGPs use
  SNR-calibrated noise, `sigma = sqrt(||X beta||^2 / (n * SNR))`.
- [trex_scr_bb_common.py](trex_scr_bb_common.py) — Biobank helpers: the control
  builder `_build_biobank_ctrl`, the module-level workers
  `_biobank_single_worker` / `_biobank_multi_worker`, the parallel runner
  `run_mc_biobank`, the per-SNR aggregator `accumulate_snr` (conditional FDR/TPR
  + unconditional Usage % and estimated FDR), and the table writer
  `save_and_print_biobank_mc`.

Both runners follow the suite convention: workers are module-level functions and
receive only picklable flat dicts; the `trex_selector_neo` enums and controls
are rebuilt inside each worker, never sent across a process boundary. Every
trial constructs its selector with `seed = -1` (hardware entropy per trial),
which the library requires for valid Monte Carlo FDR estimates.

---

## Control parameters

- `ScreenTRexControlParameter` — `trex_method`, `use_bootstrap_CI`, `R_boot`
  (Bootstrap-CI replicates, default 1000), `ci_grid_step`, and the DA knobs
  `rho_thr_DA`, `n_blocks`, `cor_coef`. It nests its own `trex_ctrl`
  (a `TRexControlParameter`); the Python `TRexScreeningSelector` also accepts a
  separate `trex_control` argument for convenience.
- `TRexControlParameter` — `solver_type`, `K`, `use_memory_mapping`,
  `solver_params.rho_afs`, `tloop_stagnation_stop`, `tloop_max_stagnant_steps`,
  … (the same object as the classical `trex/` demos).
- `BiobankScreenTRexControl` — `lower_bound_FDR`, `upper_bound_FDR` (the
  acceptable-FDR window for Algorithm 1's Ordinary → Bootstrap-CI → T-Rex-fallback
  routing), and `target_FDR_trex` (the fallback target). It nests
  `trex_screen_ctrl` (a `ScreenTRexControlParameter`, which in turn nests
  `trex_ctrl`).

`TRexBiobankScreeningSelector.select()` returns a single
`BiobankScreenTRexResult` for a 1-D response `y`, or a
`list[BiobankScreenTRexResult]` for a 2-D response matrix `Y`. Each result
exposes `phenotype_index`, `method_used` (`"Screen-TRex (ordinary)"`,
`"Screen-TRex (bootstrap-CI)"`, `"T-Rex (fallback)"`), `estimated_FDR`,
`estimated_FDR_screen_ordinary`, `estimated_FDR_screen_bootstrap`,
`selected_indices` (and `selected_indices_screen_ordinary` /
`selected_indices_screen_bootstrap`), and `used_fallback_trex`. This differs
from the R binding, which returns a list with a `$statistics` data frame; the
Python binding returns the per-phenotype result objects directly (mirroring the
C++ `BiobankScreenTRexResult` vector).

---

## Monte Carlo counts

The C++ demos use `num_MC = 200` per grid point. To keep each Python demo to
roughly a minute or two, the committed MC counts here are **downscaled**
(100 for demos 01/02, 50 for demos 03/06, 60/20 for demo 04, 40/15 for demo 05);
the biobank demos additionally use `R_boot = 500` instead of 1000. Each demo
header states its own counts. Override them without editing the sources via the
environment variables `SCR_NUM_MC`, `SCR_NUM_MC_MULTI` (biobank multi-phenotype),
and `SCR_NUM_WORKERS`; e.g. `SCR_NUM_MC=3 python <demo>.py` for a quick check.
Results are statistically comparable to the C++/R reference, **not**
bit-identical: the three use different RNG streams, so seeds are mirrored as
labels only.

---

## Running and outputs

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_01_mc_sim_screen_trex/demo_trex_scr_01_mc_sim_screen_trex.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_02_mc_sim_screen_trex_mmap/demo_trex_scr_02_mc_sim_screen_trex_mmap.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_03_mc_sim_correlated/demo_trex_scr_03_mc_sim_correlated.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_04_mc_sim_biobank/demo_trex_scr_04_mc_sim_biobank.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_05_mc_sim_biobank_mmap/demo_trex_scr_05_mc_sim_biobank_mmap.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_06_mc_sim_solvers/demo_trex_scr_06_mc_sim_solvers.py
```

Demos run from any working directory. Each writes into its own
`simulation_results/data/` subfolder:

- `.txt` — human-readable summary tables,
- `.csv` — tidy long-format tables (`method, metric, SNR|rho, value` for the
  screening demos; `method, metric, snr, value` with a `Usage` metric for the
  biobank demos) for plotting and post-processing.

Until a demo is run, `simulation_results/` holds only a `.gitkeep`.

---

**Last updated**: 2026-07-21
