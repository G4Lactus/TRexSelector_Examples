# Demo 06: DA-TRex+AKG (A Priori Known Groups) on Multi-Level Nested Data (Python)

Monte-Carlo simulations for **DA-TRex+AKG** — the dependency-aware T-Rex
selector with **a priori known groups**: the true group structure is supplied
to the selector as a hard **constraint** on its dependency tree (constrained
sub-clustering), instead of letting HAC discover the structure from scratch.
The data-generating process is a three-level nested latent-factor model with a
non-exchangeable Toeplitz leaf layer.
Common setup: n=300, p=1000, s=10, `Random` support, tFDR=0.2, K=20 random
experiments, MC=200 per grid point; solvers TLARS / TAFS / TOMP; SNR sweep
`{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`.

Python port of the C++ demo
[cpp/trex_selector_methods/trex_da/demo_trex_da_06_mc_sim_groups](../../../../cpp/trex_selector_methods/trex_da/demo_trex_da_06_mc_sim_groups/README.md)
— same DGP, control parameters, sweep axis, console summaries, and result-file
naming.

The greedy solvers use *exchangeable tie-breaking* (`exch_tie_alpha = 0.25` for
TAFS/TOMP, `0` for TLARS); TAFS additionally runs with `rho_afs = 0.3`.

## The constrained sub-clustering (PRIOR_GROUPS since 2026-07-17)

The known three-level hierarchy is passed via `da_control.prior_groups` (when
non-empty, the `method` field is ignored). The labels act as **merge
constraints**, not as deflation neighbourhoods:

1. the finest common refinement of the supplied levels (here: the groups of
   10) partitions the variables;
2. HAC runs *within* each constraint group (`hc_linkage`, correlation
   distance) — merges across group boundaries are structurally impossible;
3. the calibration rho grid is data-driven: pooled within-group dendrogram
   heights, subsampled to `hc_grid_length`, terminated by the conservative
   rho = 1 singleton anchor;
4. the BT-style occurrence deflation and 3D (T, v, rho) calibration run
   unchanged on these tight within-group neighbourhoods.

Negative group labels are rejected by the binding. Using the raw groups as
deflation neighbourhoods directly (the pre-2026-07-17 behaviour) degenerates
for group sizes > 2 — see
`TRex_Research/documentation/Prior_Groups_Deflation_Mismatch_DA_TRex.md`.

## Data-generating process

`dgp_groups_toeplitz_leaf`: each predictor is a sum of nested latent factors
plus a correlated leaf term. Three nested grouping levels over p=1000
(groups of 10 / 50 / 250, fine -> coarse) with cumulative correlation masses
`rho_levels = {0.55, 0.25, 0.10}`; within each finest group the leaf variables
follow a Toeplitz(phi=0.5) covariance, so even finest-group members are not
exchangeable. Support: s=10 actives drawn by the `Random` policy per trial,
amplitude 3.0; SNR-calibrated noise.

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_06_mc_sim_groups/demo_trex_da_06_mc_sim_groups.py
```

Monte Carlo trials run in parallel via `ProcessPoolExecutor` (up to
`min(6, os.cpu_count())` workers). Reduce `NUM_MC` to preview.

## Output files

Written to `simulation_results/`:

- `da_trex_mc_da_groups_toeplitz_leaf.txt` / `.csv`

## What to expect

See the C++ counterpart README for the committed reference results: the greedy
solvers are FDR-controlled everywhere and win on power (TOMP up to TPR ~0.99,
TAFS ~0.97), while **TLARS remains the outlier of the suite** (FDR ~0.37-0.40
from SNR=0.5 on, flat in SNR) — the LARS path distributes credit to group
siblings, producing intermediate shadow occurrences that the symmetric
nearest-partner penalty cannot separate; a structural property of the shared
BT deflation machinery under LARS-type occurrence profiles, not of the
prior-groups path.

---

**Last updated**: 2026-07-21
