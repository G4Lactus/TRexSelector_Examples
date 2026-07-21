# DA-TRex: Dependency-Aware T-Rex Demonstration Suite — Python Demos

## Overview

This folder contains the Python example programs for **DA-TRex**, the T-Rex
Selector extended with **dependency-aware (DA) correction** for
correlated/dependent design matrices. Where the classical T-Rex selector (see
[../trex/](../trex/README.md)) assumes exchangeable, weakly correlated
predictors, DA-TRex explicitly models the correlation structure among nearby
or grouped variables so that a single true signal does not get "smeared" into
false discoveries among its correlated neighbors.

The demos mirror the authoritative C++ suite in
[cpp/trex_selector_methods/trex_da/](../../../cpp/trex_selector_methods/trex_da/README.md)
one-to-one — same data-generating processes, control parameters, sweep axes,
console summaries, and result-file naming — and cover four dependency
structures and their matching DA methods: **AR(1)** (autoregressive/Toeplitz
decay), **EQUI** (equicorrelation), **NN** (banded/nearest-neighbor), and
**BT** (binary-tree/dendrogram-based hierarchical block detection via HAC
clustering), plus the **PRIOR_GROUPS** constrained sub-clustering method for a
priori known group structures.

---

## What this folder covers

- **AR(1) Toeplitz correlation** (Demo 01): the classical case DA-TRex was
  designed for — correlation decays geometrically with column distance.
- **Block-diagonal AR(1) designs**, with and without appended white-noise
  columns (Demos 03 and 02), using BT aggregation with a sweep over HAC
  linkage methods.
- **Heavy-tailed (Student-t) block designs** (Demos 04-05): testing DA-TRex
  robustness when predictors and/or noise depart from Gaussianity.
- **Multi-level nested group structures** (Demo 06): three-level hierarchical
  latent factors with a non-exchangeable Toeplitz leaf layer, driven through
  the `prior_groups` merge constraints.
- **Banded/nearest-neighbor correlation** plus a **method-mismatch stress
  test** (Demo 07: NN data, and the NN correction applied to AR(1) data).
- **Equicorrelation** (Demo 08): compound-symmetry correlation via a single
  shared latent factor.

---

## DA-TRex specific concepts

- **`tsn.DAMethod`** — the dependency-correction strategy: `AR1`, `EQUI`,
  `NN`, `BT` (plus `PRIOR_GROUPS`, engaged automatically whenever
  `da_control.prior_groups` is non-empty). Each targets the matching
  correlation structure; DA-TRex uses it to widen or narrow the "correction
  window" around each candidate variable when accumulating dummy-based voting
  evidence.
- **BT linkage method** (`LinkageMethod` from
  `trex_selector_neo.ml_methods.clustering`) — for the `BT` method, variables
  are grouped via hierarchical agglomerative clustering; several demos sweep
  **Single**, **Complete**, and **Average** linkage to study sensitivity to
  this choice.
- **Prior groups as merge constraints** (since 2026-07-17): the `PRIOR_GROUPS`
  method sub-clusters *within* the finest common refinement of the supplied
  label vectors and calibrates over a data-driven ascending rho grid anchored
  at the conservative rho = 1 singleton point; negative labels are rejected.
- **Support placement policies**: `CappedSpread` (deterministic, evenly spaced
  with a capped gap), `Random` (redrawn per trial), `OnePerBlock` (one active
  variable per block — the natural choice for block-structured DGPs). See
  `da_sim_common.make_support_*`.
- **Base T-Rex comparison rows**: Demos 01, 07, and 08 additionally run the
  *classical* T-Rex selector (no DA correction) side by side with DA-TRex on
  the same correlated data (`TREX: <solver>` rows).
- **Exchangeable tie-breaking for greedy solvers** (`exch_tie_alpha`, since
  2026-07-15): `default_solvers()` sets `exch_tie_alpha = 0.25` for TAFS/TOMP
  so that statistically indistinguishable top candidates within highly
  correlated clusters are picked uniformly at random per random experiment —
  restoring the occurrence spread the DA deflation's FDR control relies on.
  TLARS keeps `0` (path solvers spread naturally). See HISTORY.md
  (2026-07-15) in the TRexSelector repository.
- **TAFS' AFS correlation parameter**: `default_solvers()` runs TAFS with
  `rho_afs = 0.3`; TLARS and TOMP take `0` (unused).

---

## Demo descriptions

| # | Name | DA Method | DGP | Key parameters |
|---|------|-----------|-----|-----------------|
| **01** | [AR(1)](demo_trex_da_01_mc_sim_ar1/) | `AR1` | `dgp_ar1_snr` | n=300, p=1000, s=10, rho=0.7; SNR & rho sweeps + 2D gap x rho |
| **02** | [BT + white noise](demo_trex_da_02_mc_sim_ar1_blocks_plus_white/) | `BT` | `dgp_ar1_block_white_snr` | n=300, p_total=1000, M=5, Q=5 (p_ar=25), rho=0.7; linkage x {SNR, rho, Q, M} sweeps |
| **03** | [BT](demo_trex_da_03_mc_sim_ar1_blocks/) | `BT` | `dgp_ar1_block_snr` | n=150, M=5, Q=5 (p=25, s=5), rho=0.7; linkage x {SNR, rho, Q, M} sweeps |
| **04** | [BT, heavy-tailed](demo_trex_da_04_mc_sim_ht_blocks/) | `BT` | `dgp_block_toeplitz_hvt_snr` | n=150, M=5, Q=5, rho=0.8, nu=3; Gauss/Heavy x linkage x {SNR, rho, Q, M, tFDR} |
| **05** | [BT, heavy-tailed + white](demo_trex_da_05_mc_sim_ht_blocks_plus_white/) | `BT` | `dgp_ht_block_white_snr` | n=150, p_total=500, M=5, Q=5, rho=0.8, nu=3; Gauss/Heavy x linkage x {SNR, rho, Q, M, tFDR} |
| **06** | [Prior groups (AKG)](demo_trex_da_06_mc_sim_groups/) | `PRIOR_GROUPS` (a priori known groups) | `dgp_groups_toeplitz_leaf` | n=300, p=1000, s=10; 3-level groups {10, 50, 250}, rho={0.55, 0.25, 0.10}, phi=0.5 |
| **07** | [NN + mismatch test](demo_trex_da_07_mc_sim_nn/) | `NN` | `dgp_nn_snr`, `dgp_ar1_snr` | n=300, p=1000, kappa=3, rho=0.7; SNR sweep + 2D kappa x rho + 2D SNR x rho on AR(1) data (NN misspecified) |
| **08** | [Equi](demo_trex_da_08_mc_sim_equi_and_bt/) | `EQUI` | `dgp_equi_snr` | n=300, p=1000, s=10; rho sweep x SNR in {0.5, 1, 2, 5} |

---

## Folder contents

```txt
trex_da/
  ├── README.md
  ├── dgp_generators.py       ┐  suite-level shared modules
  ├── da_sim_common.py        ┘  (imported by every demo)
  ├── demo_trex_da_01_mc_sim_ar1/
  ├── demo_trex_da_02_mc_sim_ar1_blocks_plus_white/
  ├── demo_trex_da_03_mc_sim_ar1_blocks/
  ├── demo_trex_da_04_mc_sim_ht_blocks/
  ├── demo_trex_da_05_mc_sim_ht_blocks_plus_white/
  ├── demo_trex_da_06_mc_sim_groups/
  ├── demo_trex_da_07_mc_sim_nn/
  └── demo_trex_da_08_mc_sim_equi_and_bt/
```

Each demo subfolder contains its `.py` script (named after the folder), a
`README.md`, and a `simulation_results/` folder that receives the `.txt`
(pretty-printed tables) and `.csv` (tidy long format) files written by the
demo, following the cpp naming pattern `da_trex_mc_{scenario_tag}.{txt,csv}`.
(The cpp suite splits `simulation_results/` into `data/` + `plots/`; the
Python suite has no plotting layer, so results live directly in
`simulation_results/`.)

Shared infrastructure:

- [dgp_generators.py](dgp_generators.py) — the DA data-generating processes
  (AR(1); equicorrelated; NN/MA(kappa); AR(1)-block, with and without a white
  block; heavy-tailed blocks, with and without a heavy-tailed white block;
  prior-groups with Toeplitz leaves), all SNR-calibrated with 0-based support,
  mirroring `cpp/.../dgp_generators.hpp` one generator per cpp DGP.
- [da_sim_common.py](da_sim_common.py) — support generators, the default
  solver set (TLARS/TAFS/TOMP with per-solver `exch_tie_alpha` / `rho_afs` /
  stagnation), parallel MC runners for DA-TRex and base T-Rex, generic
  parameter sweeps, and the console/TXT/CSV writers, mirroring
  `cpp/.../trex_da_sim_common.hpp`.

---

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_01_mc_sim_ar1/demo_trex_da_01_mc_sim_ar1.py
```

Each demo exposes `RUN_PART_*` flags at the top; all parts are enabled by
default (mirroring the cpp `if (true)` toggles) and each runs `NUM_MC = 200`
Monte Carlo trials per grid point — reduce `NUM_MC` to preview. Each demo
inserts both its own directory and its parent suite directory onto `sys.path`,
so the shared `dgp_generators` / `da_sim_common` modules resolve from the
nested location and are importable by the spawned workers — so demos must be
run as files, not piped via `python -`. Monte Carlo trials run in parallel via
`concurrent.futures.ProcessPoolExecutor` (one persistent pool per demo run,
up to `min(6, os.cpu_count())` workers).

Seed spaces mirror the cpp suite: DGP seed = base + grid offset + trial,
random support seed = trial seed + 500000, selector seed = -1 (hardware
entropy, required for valid MC FDR estimates).

Indices are **0-based** (Python convention); the R counterparts are 1-based.

---

## Start here

1. **Demo 01 — AR(1)**: the foundational scenario, including a direct
   DA-TRex vs. base-T-Rex comparison and a "gap x rho" sweep showing when
   correlated neighbors start hurting the classical selector.
2. **Demo 03 — BT on block AR(1)**: introduces the `BT` method and HAC
   linkage sweep on a clean block-diagonal design.
3. **Demo 06 — Prior groups**: the most structurally rich DGP (three nested
   levels of correlation), driven through the constrained sub-clustering.

Cross-check against the C++ suite (committed reference data and figures) in
[cpp/trex_selector_methods/trex_da/](../../../cpp/trex_selector_methods/trex_da/README.md)
and the R suite in `R/trex_selector_methods/trex_da/`. Note that the R suite
keeps its own legacy numbering (e.g. the NN demos are `03`/`03b` there), so
match demos by name, not number.

---

**Last updated**: 2026-07-21
