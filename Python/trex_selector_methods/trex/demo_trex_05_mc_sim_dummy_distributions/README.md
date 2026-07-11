# Demo 05: Dummy Distribution Comparison (Python)

## Purpose

Compare the T-Rex **dummy distributions** (the distribution from which the
injected dummy/null variables are drawn) across three base solvers — TLARS
(equiangular LARS path), TOMP, and TAFS with `rho_afs = 1.0` (both greedy) —
with the L-loop strategy fixed at STANDARD. Since dummies are centered and
L2-normalized before entering the solver, their *scale* is immaterial; the demo
probes whether the *shape* (tail weight, discreteness, skewness, sparsity,
cross-column dependence) affects FDR, TPR, calibrated dummy multiplier `L`, and
stopping time `T`. Python port of
`cpp/trex_selector_methods/trex/demo_trex_05_mc_sim_dummy_distributions`.

## DGP / Data

- `n = 300`, `p = 1000` (high-dimensional), `s = 10`
- Random support redrawn per trial by default (`block_support=False`); a
  contiguous block support `{0, ..., s-1}` is available via `block_support=True`
- Fixed coefficients `beta_j = 1` (`rnd_coef=False`)
- SNR grid: `{0.1, 0.2, ..., 2.0}` plus `5.0` (21 values)
- `num_MC = 10` trials per solver x distribution x SNR (`_NUM_MC = 10`, mirrors C++)
- Base solvers (`_SOLVERS`, outer sweep): TLARS, TOMP, TAFS (`rho_afs = 1.0`);
  L-loop strategy: STANDARD (fixed)

## Distributions

`DISTRIBUTIONS` defines **12 rows** as plain picklable spec dicts, resolved to
`trex_selector_neo.DummyDistribution` objects inside each worker by
`trex_sim_common._make_dummy_distribution()`:

Normal, Uniform, Rademacher, StudentT_df3, StudentT_df5, Laplace, Gumbel,
Triangle, Logistic, Mammen, SparseRad_s0.1 (constrained sparse Rademacher,
`s = 0.1`), and UnifSphere_d5 (uniform on the unit 5-sphere; consecutive
5-column groups are dependent — a deliberate probe of the exchangeability
axis).

**Note on UniformSphere**: the library requests dummies in multiples of `p`, so
the sphere dimension must divide `p`; the default `dim = 3` would throw at
`p = 1000`, hence `dim = 5`.

## What It Computes

`trx_05_dummy_distributions()` sweeps each solver x distribution x SNR via
`run_mc_trex` with a per-distribution control (`K = 20`, `max_dummy_multiplier = 10`,
`lloop_strategy = "STANDARD"`, `dummy_distribution` from the spec). With
`track_L_T=True` it records mean FDR, TPR, Avg L, and Avg T per distribution x
SNR.

## Expected pattern

The C++ reference run shows all i.i.d. unit-variance rows FDR-equivalent at or
below `tFDR = 0.1` (isolated exceedances at very low SNR are small-`num_MC`
noise), statistically indistinguishable TPR, and no systematic deviation even
for the dependent UnifSphere_d5 rows.

## Imports

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `trex_sim_common` module. The package imports as
`trex_selector_neo`.

## Output

Writes results to this demo's own `simulation_results/data/` folder via
`save_mc_results`, with one file pair per base solver, stems
`demo_trex_05_dummy_distributions_results_n300_p1000_{tlars,tomp,tafs}_random_support`
(aligned `.txt` table plus tidy `.csv`; the "solver" column holds the
distribution name — the base solver is encoded in the file name).

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex/demo_trex_05_mc_sim_dummy_distributions/demo_trex_05_mc_sim_dummy_distributions.py
```

The MC loop is parallelized with Python `multiprocessing` (spawn start method)
using `_NUM_WORKERS = 6`.

**Last updated**: 2026-07-11
