# Demo 05: Dummy Distribution Comparison

## Purpose

Compare the T-Rex **dummy distributions** — the distribution from which the injected
dummy (null) variables are drawn — across three base solvers
(**TLARS** — equiangular LARS path; **TOMP** and **TAFS** with `rho_afs = 1.0` —
greedy forward selection), with the L-loop strategy fixed at **STANDARD**. Since dummies are centered and
L2-normalized before entering the solver, their *scale* is immaterial; the demo probes
whether the *shape* (tail weight, discreteness, skewness, sparsity, cross-column
dependence) affects FDR, TPR, and the calibrated dummy multiplier `L` / stopping time `T`.

## Data-generating process

`dgp_gauss_snr` — i.i.d. Gaussian design, `y = X beta + eps` calibrated to SNR.
High-dimensional config: `n = 300`, `p = 1000`, `s = 10`, `tFDR = 0.1`, fixed
coefficients `beta_j = 1`. The support is a random size-`s` subset redrawn per trial
via `make_support_random(s, p, seed = trial_seed)` (the `block_support = TRUE`
contiguous-block variant is available but its call is commented out). SNR grid
`seq(0.1, 2.0, by = 0.1)` plus `5.0` (21 values), `num_MC = 10` per solver x distribution x SNR
(mirrors the C++ demo).

The solver sweep (`SOLVERS`) is the outer loop; one `.txt` + `.csv` file pair is
written per solver, with the lowercased solver name in the stem
(`..._{tlars,tomp,tafs}_random_support`).

## Distributions compared

`DISTRIBUTIONS` — **12 rows**, each built with `trex_dummy_distribution()` and passed
to `trex_control(dummy_distribution = ...)`:

- Normal — N(0, 1); the reference choice.
- Uniform — U(-sqrt(3), sqrt(3)); bounded, light tails.
- Rademacher — {-1, +1} equiprobable; discrete two-point.
- StudentT_df3 / StudentT_df5 — heavy tails (unit-variance scaled).
- Laplace — double exponential; heavier-than-normal tails.
- Gumbel — extreme-value; skewed (centered at 0).
- Triangle — bounded, symmetric; (-sqrt(6), 0, sqrt(6)).
- Logistic — tail weight between normal and Laplace.
- Mammen — asymmetric golden-ratio two-point distribution.
- SparseRad_s0.1 — constrained sparse Rademacher, sparsity `s = 0.1`.
- UnifSphere_d5 — uniform on the unit 5-sphere; consecutive 5-column groups are
  dependent (norm constraint), a deliberate probe of the exchangeability axis.

**Note on UniformSphere**: the library requests dummies in multiples of `p`, so the
sphere dimension must divide `p`; the default `dim = 3` would throw at `p = 1000`,
hence `dim = 5`.

Shared control per solver x distribution: `K = 20`,
`max_dummy_multiplier = 10`, `use_max_T_stop = TRUE`, `lloop_strategy = "STANDARD"`,
`dummy_distribution` swept. Reports averaged FDR, TPR, Avg L, Avg T per distribution x
SNR (`track_L_T = TRUE`). MC trials run in parallel over a hardcoded `num_cores <- 6L`.

## Expected pattern

The C++ reference run shows all i.i.d. unit-variance rows FDR-equivalent at or below
`tFDR` (isolated exceedances at very low SNR are small-`num_MC` noise), statistically
indistinguishable TPR, and no systematic deviation even for the dependent
UnifSphere_d5 rows.

## Running

```bash
Rscript R/trex_selector_methods/trex/demo_trex_05_mc_sim_dummy_distributions/demo_trex_05_mc_sim_dummy_distributions.R
```

Results are written to this demo's own `simulation_results/data/` folder as a `.txt` table
and a tidy `.csv` (the "solver" column holds the distribution name), stems
`demo_trex_05_dummy_distributions_results_n300_p1000_{tlars,tomp,tafs}_random_support`. The script sources
`../../support_generators.R`, `../../simulation_utils.R`, and `../trex_sim_utils.R`.

**Last updated**: 2026-07-11
