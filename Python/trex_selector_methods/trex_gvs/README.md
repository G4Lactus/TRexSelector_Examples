# T-Rex+GVS: Grouped Variable Selection — Python Demos

## Overview

This folder contains the Python ports of the **T-Rex Selector with Grouped Variable Selection (T-Rex+GVS)**
demonstration suite according to [[1,2]](#references), using `trex_selector_neo`'s `TRexGVSSelector`.
They mirror the authoritative C++ suite in
[cpp/trex_selector_methods/trex_gvs/](../../../cpp/trex_selector_methods/trex_gvs/README.md) one-to-one
(folder names, scenarios, sweep grids, console tables, and result-file names all match).

T-Rex+GVS extends the classical T-Rex selector (see [../trex/](../trex/README.md)) to settings where predictors
 form known or discoverable **groups** of correlated variables.
 Examples are genes in the same pathway, or sensor channels driven by a shared latent signal.
 The selector aggregates the evidence of the random experiments at the group level while it runs, so that correlated
 variables reinforce each other instead of splitting the vote.
 However, the selection is based on individual variables, and that is also the level at which FDR and TPR are
 evaluated (see [What is actually measured](#what-is-actually-measured-in-these-demos)).

The demos in this folder are designed to help users understand:

1. how grouped-signal recovery behaves under different within-group correlation structures (equicorrelated, AR(1),
    ARMA, heavy-tailed),
2. how group *discovery* via hierarchical agglomerative clustering (HAC) compares to using known ("oracle") groups,
3. how the elastic-net-based (`EN`) and informed-elastic-net-based (`IEN`) group solvers behave under scattered
    layouts, correlation traps, and negative correlation structure.

If you are new to this folder, start with **Demo 01** for the canonical equicorrelated-blocks scenario.

---

## What this folder covers

All demos study sparse variable selection in a Gaussian (or near-Gaussian) linear model where the active variables
 cluster into groups:

- **equicorrelated blocks** (Demos 01–04, 08): variables within a group share a common latent factor,
- **scattered groups** (Demo 02): group members are permuted across the column space instead of being contiguous,
- **gapped, randomly-ordered blocks with an inactive trap** (Demo 03, and reused by Demos 05–07): four unequal-size
   blocks, three active and one inactive equicorrelated trap,
- **negative-correlation structure** (Demo 04): groups with sign-flipped halves and spatially separated inactive traps,
- **non-Gaussian / non-equicorrelated structure** (Demos 05–07): Student-t(3) heavy tails, AR(1) Toeplitz covariance,
   and heterogeneous ARMA processes per block,
- **value of prior group information** (Demo 08): a larger 100-block benchmark comparing HAC-discovered groups
   against known ("oracle") group labels, sweeping the within-block correlation across the HAC discoverability
   boundary.

| # | Demo | Scenario | Output stem | cpp counterpart |
|---|---|---|---|---|
| 01 | [demo_trex_gvs_01_mc_sim_hastie_en_blocks](demo_trex_gvs_01_mc_sim_hastie_en_blocks/) | Hastie: 3 equicorrelated active groups of 50 (s = 150); 2-D grid + lambda2/linkage studies | `gvs_Hastie_*` | [cpp demo 01](../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_01_mc_sim_hastie_en_blocks/README.md) |
| 02 | [demo_trex_gvs_02_mc_sim_scattered_blocks](demo_trex_gvs_02_mc_sim_scattered_blocks/) | Same 3 groups, members permuted across columns | `gvs_scattered_grouped_2d` | [cpp demo 02](../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_02_mc_sim_scattered_blocks/README.md) |
| 03 | [demo_trex_gvs_03_mc_sim_mixed_blocks](demo_trex_gvs_03_mc_sim_mixed_blocks/) | 3 active blocks (20/50/80) + 1 trap (65), random order, noise gaps | `gvs_mixed_blocks_2d` | [cpp demo 03](../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_03_mc_sim_mixed_blocks/README.md) |
| 04 | [demo_trex_gvs_04_mc_sim_neg_traps](demo_trex_gvs_04_mc_sim_neg_traps/) | Active +Z/−Z group + two inactive +Z/−Z traps (s = 100) | `gvs_Neg-Traps_2d` | [cpp demo 04](../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_04_mc_sim_neg_traps/README.md) |
| 05 | [demo_trex_gvs_05_mc_sim_t3_blocks](demo_trex_gvs_05_mc_sim_t3_blocks/) | Demo-03 geometry with Student-t(3) factors/noise | `gvs_t3_blocks_2d` | [cpp demo 05](../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_05_mc_sim_t3_blocks/README.md) |
| 06 | [demo_trex_gvs_06_mc_sim_ar1_blocks](demo_trex_gvs_06_mc_sim_ar1_blocks/) | Demo-03 geometry with within-block Cor = rho^\|j−k\| | `gvs_ar1_blocks_2d` | [cpp demo 06](../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_06_mc_sim_ar1_blocks/README.md) |
| 07 | [demo_trex_gvs_07_mc_sim_arma_blocks](demo_trex_gvs_07_mc_sim_arma_blocks/) | AR(1) / MA(3) / ARMA(2,1) active blocks + AR(1) trap; ar_coef axis | `gvs_arma_blocks_2d` | [cpp demo 07](../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_07_mc_sim_arma_blocks/README.md) |
| 08 | [demo_trex_gvs_08_mc_sim_block_bench](demo_trex_gvs_08_mc_sim_block_bench/) | 100 blocks × 5, 10 active; HAC vs oracle groups (M1–M4), WHOLE + INDIVIDUAL truths | `gvs_block_bench_*` | [cpp demo 08](../../../cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_08_mc_sim_block_bench/README.md) |

Demos 01–07 compare the three solver variants **TENET** (EN), **TENET_AUG** (EN, row-augmented solver), and
**TIENET_AUG** (IEN) on a 2-D SNR × correlation grid (fixed settings: `n = 200`, `p = 500`, `tFDR = 0.1`,
`K = 20`, `corr_max = 0.98`, 200 MC trials per grid point). Demo 01 adds `lambda2_method` and `hc_linkage`
robustness checks; Demo 08 compares EN/IEN × HAC/oracle (M1–M4) at `corr_max = 0.5`.

---

## What is actually measured in these demos

> **Selection and evaluation are at the level of individual variables.** Although T-Rex+GVS
> *aggregates the evidence of the random experiments at the group level* while selecting, what the selector returns
> is a set of individual variable indices, and every FDR/TPR number reported in Demos 01–07 is
> **coordinate-level**: a false discovery is a selected variable that is not truly active.
>
> A genuinely **group-level FDR** — treating a whole group as the unit of discovery — is a
> different quantity and needs its own definition of what counts as a false discovery, together with a matching
> notion of the group null hypothesis; see [[3]](#references) for such a formulation in the knockoff setting. It is
> *not* obtainable by reinterpreting the coordinate-level numbers.

Demo 08 is the only demo that additionally reports block-level diagnostics (block hit rate, block FDR, full-block
 rate, block purity). These are descriptive only: the selector controls the coordinate-level FDR, and no guarantee
 attaches to any block-level rate.

For the statistical model, the within-group correlation models, and the FDR/TPR definitions, see the
[C++ suite README](../../../cpp/trex_selector_methods/trex_gvs/README.md#statistical-model-and-targets) — the
Python DGPs implement the same constructions (with NumPy's PRNG, so per-cell numbers differ slightly from the
committed C++ results while reproducing the same qualitative patterns).

---

## T-Rex+GVS specific concepts

Grouped selection introduces a few additional control parameters and types beyond the classical T-Rex selector
 (all accessed via `import trex_selector_neo as tsn`, e.g. `tsn.TRexGVSControlParameter`):

- **`tsn.GVSType`** — the group-level base method:
  - `EN`: elastic-net-based grouped selection [[1]](#references),
  - `IEN`: informed (iterative) elastic-net-based grouped selection [[2]](#references).
- **`tsn.ENSolverType`** (only applies to `GVSType.EN`) — which elastic-net solver drives the random experiments:
  - `TENET`: the base Gram-based elastic-net solver,
  - `TENET_AUG`: the augmented-LASSO formulation of the elastic net.
- **`tsn.LambdaSelectionMethod`** — how the elastic-net penalty $\lambda_2$ is chosen: `CV_1SE_CCD`
   (cross-validation, 1-SE rule, via coordinate descent) or `CV_1SE_SVD` (via SVD).
- **`prior_groups`** — the 0-based per-column group-ID vector. When groups are known in advance (an "oracle" run),
   this vector is taken directly from the data-generating process. When groups must be discovered from the data,
   HAC clustering (with a `corr_max` correlation threshold and a chosen linkage method) is used to build
   `prior_groups` automatically — this is what Demo 08 compares directly against oracle groups.

The `TRexGVSSelector` wrapper derives the required solver from `gvs_type` / `en_solver`
(EN → TENET/TENET_AUG, IEN → TIENET_AUG), so the demos never set `trex_control.solver_type` themselves.

---

## Folder contents

```txt
trex_gvs/
  ├── README.md
  ├── trex_gvs_dgps.py                # all scenario DGPs (mirrors trex_gvs_dgps.hpp)
  ├── gvs_sim_common.py               # selector adapter, parallel MC runners,
  │                                   # table printers, .txt/.csv writers
  │                                   # (mirrors trex_gvs_simulation_utils.hpp)
  │
  ├── demo_trex_gvs_01_mc_sim_hastie_en_blocks/
  │   ├── demo_trex_gvs_01_mc_sim_hastie_en_blocks.py   # the demo script
  │   ├── README.md                    # scenario description and interpretation
  │   └── simulation_results/
  │       └── data/                    # .txt summaries + .csv tables
  │
  │   ... demos 02–08 follow the same layout ...
  │
  └── demo_trex_gvs_08_mc_sim_block_bench/
      ├── demo_trex_gvs_08_mc_sim_block_bench.py
      ├── README.md
      └── simulation_results/
          └── data/
```

Shared infrastructure:

- [trex_gvs_dgps.py](trex_gvs_dgps.py) — all scenario data-generating processes in one module:
  `dgp_hastie_snr`, `dgp_scattered_grouped_snr`, `dgp_mixed_blocks_snr`, `dgp_neg_traps_snr`,
  `dgp_t3_equi_blocks_snr`, `dgp_ar1_blocks_snr`, `dgp_arma_blocks_snr`, and `dgp_block_equicorr`.
  All are SNR-calibrated (`sigma_y = sqrt(var(signal) / snr)`, `ddof=1`; the block benchmark uses
  the biased `ddof=0` denominator, matching the cpp DGP). Support / group indices are 0-based.
- [gvs_sim_common.py](gvs_sim_common.py) — the `build_gvs_selector` adapter (EN / EN_AUG / IEN, HAC or
  oracle `prior_groups`), the parallel MC runners `run_gvs_mc_trials` / `run_gvs_2d_sweep` /
  `run_block_mc` (`ProcessPoolExecutor`, mirroring the cpp OpenMP loops), and the table printers and
  `.txt`/`.csv` writers whose output matches the cpp layout and schemas.

---

## Running a demo

Each demo is run directly with the Python interpreter (with `trex_selector_neo`, `numpy`, and
`pandas` installed — see the [Python examples README](../../README.md)):

```bash
cd Python/trex_selector_methods/trex_gvs/demo_trex_gvs_01_mc_sim_hastie_en_blocks
python demo_trex_gvs_01_mc_sim_hastie_en_blocks.py       # full run (hours)
python demo_trex_gvs_01_mc_sim_hastie_en_blocks.py 8     # with 8 worker processes
```

Notes:

- Every demo writes a human-readable `.txt` summary and a tidy `.csv` table into its local
  `simulation_results/data/` folder, using the same file names and CSV schemas as the C++ suite
  (so the cpp-side `trex_gvs_plt_utils.py` can render the same figures from the Python CSVs).
- Each demo inserts its own directory and the parent suite directory onto `sys.path`, so the shared
  `trex_gvs_dgps` / `gvs_sim_common` modules resolve from the nested location and are re-importable
  by the `spawn`-launched Monte Carlo workers — run the demos **as files**, not piped via `python -`.
- Monte Carlo trials are parallelized with `concurrent.futures.ProcessPoolExecutor`. The optional
  first CLI argument sets the number of worker processes (default `min(6, os.cpu_count())`).
- The committed settings match the C++ suite (200 MC trials per grid point) and take hours per demo;
  for a quick smoke test, temporarily reduce `num_MC` and the sweep grids at the top of a script.
- Seeding mirrors the cpp MC loops: trial i uses data seed `base_seed + i`, which is also passed as
  the selector seed (and as `cv_seed` in demos 01–07); demo 08 staggers seeds per grid cell.
- Indices are **0-based** (Python convention); the R counterparts are 1-based.

---

## References

1. Machkour, J., Muma, M., & Palomar, D. P., "False Discovery Rate Control for Grouped Variable Selection
   in High-Dimensional Linear Models using the T-Knock Filter.", European Signal Processing Conference (EUSIPCO), 2022,
    pp. 892–896, EURASIP.
    [DOI-Link](https://doi.org/10.23919/EUSIPCO55093.2022.9909883)
2. Machkour, J., Muma, M., & Palomar, D. P., "The Informed Elastic Net for Fast Grouped Variable Selection and
   FDR Control in Genomics Research.", Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP),
    2023, pp. 466–470, IEEE.
    [DOI-Link](https://doi.org/10.1109/CAMSAP58249.2023.10403489)
3. Dai, R., & Barber, R. F., "The knockoff filter for FDR control in group-sparse and multitask regression."
   Proceedings of the 33rd International Conference on Machine Learning (ICML), vol. 48, pp. 1851–1859, 2016.
   [DOI-Link](http://proceedings.mlr.press/v48/daia16.pdf)

---

**Last updated**: 2026-07-21
