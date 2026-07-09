# Screen-TRex Selector: Python Demonstration Suite

## Purpose

Python example programs for **Screen-TRex**, the T-Rex Selector extended with a
variable-**screening** step for ultra-high-dimensional data. Screen-TRex is
built for the regime where `p` is so large that running the classical T-Rex
selector directly is impractical (genome- or biobank-scale designs): it first
thresholds a dummy-based voting statistic `Phi_j` to screen down to a small
candidate set, using either the **Ordinary** majority-vote rule
(`Phi_j > 0.5`) or a **Bootstrap-CI** rule, plus **dependency-aware (DA)**
variants that recover power under correlated designs.

```python
import trex_selector_neo as ts

# Ordinary / DA screening
sc = ts.ScreenTRexControlParameter()
sc.trex_method = ts.ScreenTRexMethod.TREX          # or a TREX_DA_* variant
tc = ts.TRexControlParameter(); tc.K = 20
sel = ts.TRexScreeningSelector(X, y, screen_control=sc, trex_control=tc,
                               seed=42, verbose=False)
res = sel.select()
sel.selected_indices          # 0-based screened support (numpy array)
res.estimated_FDR             # procedure's self-reported FDR
res.estimated_correlation     # estimated rho (DA variants)

from trex_selector_neo.utils import compute_fdp, compute_tpp
compute_fdp(sel.selected_indices.tolist(), true_support)
compute_tpp(sel.selected_indices.tolist(), true_support)

# Biobank / Algorithm 1 (many phenotypes, one shared X)
bc = ts.BiobankScreenTRexControl()
bc.lower_bound_FDR = 0.05; bc.upper_bound_FDR = 0.15; bc.target_FDR_trex = 0.10
bio = ts.TRexBiobankScreeningSelector(X, Y, bio_ctrl=bc, seed=42, verbose=False)
results = bio.select()        # BiobankScreenTRexResult, or a list for m-D Y
results[0].method_used        # routing decision per phenotype
results[0].selected_indices   # per-phenotype 0-based selected indices
```

The suite mirrors the C++ demos in
[cpp/trex_selector_methods/trex_screening/](../../../cpp/trex_selector_methods/trex_screening/README.md)
(see that README for the statistical model, notation, `Phi`, DA screening, and
Algorithm 1) and the R demos in
[R/trex_selector_methods/trex_screening/](../../../R/trex_selector_methods/trex_screening/README.md).
Support sets and selected indices are **0-based** in Python (matching the C++
sources; the R suite is 1-based).

The `ScreenTRexMethod` values accepted by `ScreenTRexControlParameter.trex_method`
are `ts.ScreenTRexMethod.TREX` (ordinary voting-based screening), `.TREX_DA_AR1`,
`.TREX_DA_EQUI`, and `.TREX_DA_BLOCK_EQUI` (dependency-aware variants).

The folder layout mirrors the C++ suite one-to-one: one subfolder per demo,
plus the two shared helper modules at the suite level.

```txt
trex_screening/
  ├── README.md
  ├── trex_scr_common.py              <- DGPs, control factories, MC runner, table output
  ├── trex_scr_bb_common.py           <- Biobank (Algorithm 1) helpers + MC table output
  ├── demo_trex_scr_01_screen/
  │   ├── demo_trex_scr_01_screen.py
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
| 01 | [demo_trex_scr_01_screen/](demo_trex_scr_01_screen/) | In-memory Ordinary vs. Bootstrap-CI screening: two single runs (fixed support, SNR 5) plus an 8-point SNR-sweep Monte Carlo reporting FDR / TPR / Estimated FDR | [demo_trex_scr_01_screen_trex/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_01_screen_trex/) | `TRexScreeningSelector`, `ScreenTRexControlParameter.use_bootstrap_CI` |
| 02 | [demo_trex_scr_02_screen_mmap/](demo_trex_scr_02_screen_mmap/) | Same screening scenario with the memory-mapped pipeline: Part A in-memory X + `use_memory_mapping=True` (D-mmap); Part B fully memory-mapped X via `numpy_to_memmap()`; Part C mmap MC SNR sweep | [demo_trex_scr_02_screen_trex_mmap/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_02_screen_trex_mmap/) | `numpy_to_memmap`, `TRexControlParameter.use_memory_mapping` |
| 03 | [demo_trex_scr_03_screen_correlated/](demo_trex_scr_03_screen_correlated/) | Ordinary vs. DA screening on correlated designs: 4 single runs (Ordinary/DA-AR1/DA-EQUI/DA-BLOCK-EQUI) + 5 Monte Carlo sweeps (AR(1)/equi SNR sweeps and AR(1)/equi/block-equi rho sweeps) | [demo_trex_scr_03_screen_trex_correlated/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_03_screen_trex_correlated/) | `ScreenTRexMethod.TREX_DA_*`, `estimated_correlation` |
| 04 | [demo_trex_scr_04_biobank_inmem/](demo_trex_scr_04_biobank_inmem/) | Biobank "Algorithm 1" (in-memory): single phenotype, 5-phenotype routing table + summary, and MC SNR sweeps (single and multi-phenotype) reporting per-method Usage % / FDR / TPR / Est. FDR | [demo_trex_scr_04_biobank_screen_trex_inmem/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_04_biobank_screen_trex_inmem/) | `TRexBiobankScreeningSelector`, `BiobankScreenTRexControl` |
| 05 | [demo_trex_scr_05_biobank_mmap/](demo_trex_scr_05_biobank_mmap/) | Same biobank workflow with the memory-mapped pipeline: Part A in-memory X + D-mmap; Part B single + multi phenotype on a fully memory-mapped X; Part C mmap MC SNR sweep | [demo_trex_scr_05_biobank_screen_trex_mmap/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_05_biobank_screen_trex_mmap/) | `TRexBiobankScreeningSelector` + mmap X |
| 06 | [demo_trex_scr_06_screen_solvers/](demo_trex_scr_06_screen_solvers/) | Solver-backend comparison under screening: 8-point SNR sweep for TLARS / TAFS(rho_afs=0.3) / TOMP, Part 1 Ordinary only (3 series), Part 2 x {Ordinary, Bootstrap-CI} (6 series) | [demo_trex_scr_06_screen_trex_solvers/](../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_trex_solvers/) | `TRexControlParameter.solver_type`, `solver_params.rho_afs` |

---

## Shared modules

Both live at the suite level and mirror the two C++ demo-internal headers
(`demo_trex_scr_common.hpp`, `demo_trex_scr_bb_common.hpp`). Each demo resolves
them via a `sys.path` bootstrap that inserts both its own directory and the
parent (suite) directory, so the shared modules import cleanly from the nested
demo location.

- [trex_scr_common.py](trex_scr_common.py) — DGPs (`dgp_iid_snr`, `dgp_ar1`,
  `dgp_equi`, `dgp_block_equi`, `make_beta`), control factories
  (`make_scr_trex_ctrl`, `make_screen_ctrl`, `default_scr_methods`), the
  single-run printer `print_scr_result`, the Monte Carlo runner
  `run_mc_screen`, and the table writer `save_and_print_scr_mc`
  (console + `.txt` + tidy `.csv`). All DGPs use SNR-calibrated noise,
  `sigma = sqrt(||X beta||^2 / (n * SNR))`.
- [trex_scr_bb_common.py](trex_scr_bb_common.py) — Biobank helpers:
  `make_biobank_ctrl`, `print_biobank_single`, `print_biobank_table`,
  `print_biobank_summary`, and `save_and_print_biobank_mc` (with a per-method
  Usage (%) row).

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

The C++ demos use `num_MC = 500`. To keep each Python demo to roughly a minute
or two single-process, the MC counts here are reduced (60 for demo 01, 40 for
demo 02, 10 for demo 03's five sweeps, 40/15 for demo 04, 30 for demo 05, 20 for
demo 06) — each demo header states its own count. Results are statistically
comparable to the C++/R reference, **not** bit-identical: the three use
different RNG streams, so seeds are mirrored as labels only.

---

## Running and outputs

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_01_screen/demo_trex_scr_01_screen.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_02_screen_mmap/demo_trex_scr_02_screen_mmap.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_03_screen_correlated/demo_trex_scr_03_screen_correlated.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_04_biobank_inmem/demo_trex_scr_04_biobank_inmem.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_05_biobank_mmap/demo_trex_scr_05_biobank_mmap.py
python Python/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_solvers/demo_trex_scr_06_screen_solvers.py
```

Demos run from any working directory. The Monte Carlo parts write into each
demo's own `simulation_results/` subfolder:

- `.txt` — human-readable summary tables,
- `.csv` — tidy long-format tables (`method, metric, SNR|rho, value`) for
  plotting and post-processing.

The single-run parts (and the biobank single/multi-phenotype tables) print to
the console only. The `RUN_PART_*` flags at the top of each demo toggle
individual parts.

---

**Last updated**: 2026-07-08
