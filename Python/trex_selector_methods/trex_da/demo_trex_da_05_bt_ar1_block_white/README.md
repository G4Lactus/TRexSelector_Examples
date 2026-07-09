# Demo 05: T-Rex+DA+BT on Block AR(1) + White-Noise Data (Python)

## Purpose

Extend Demo 04's block-diagonal AR(1) design by appending a large block of i.i.d.
white-noise columns, testing whether `BT`-based DA-TRex can still isolate the true
AR(1) blocks when they are a small fraction of a much larger, otherwise-uncorrelated
design. Sensitivity to the HAC linkage method is studied across four sweep
dimensions.

Python port of `R/trex_selector_methods/trex_da/demo_trex_da_05_bt_ar1_block_white`
and `cpp/.../demo_trex_da_05_mc_sim_bt_ar1_block_sweeps`.

## Data-generating process

`dgp_ar1_block_white_snr` builds M AR(1) blocks of size Q (`p_ar = M * Q`,
structured) plus `p_white = p_total - p_ar` i.i.d. N(0,1) white-noise columns, with
`p_total` fixed at 1000. Active variables (`s = M`, one per AR block) lie only in the
AR part. Base config: `n=300, p_total=1000, M=5, Q=5` (p_ar=25, p_white=975, s=5),
`rho=0.7, snr=2.0, amplitude=1.0, tFDR=0.2, K=20, num_MC=200, seed=2026`.

## DA method

`build_da_selector("BT", ..., hc_linkage=lnk)` constructs a
`trex_selector_neo.TRexDASelector` with a `TRexDAControlParameter`
(`method = DAMethod.BT`, `hc_linkage` from `LinkageMethod.{Single,Complete,Average}`)
and a `TRexControlParameter` (`K = 20`, `solver_type = SolverTypeForTRex.TLARS`). Every
sweep loops the linkage over `{single, complete, average}`.

## What it runs

All parts active by default (`RUN_PART_1..4 = True`). Each is heavy — reduce `num_MC`
to preview.

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- Part 2: rho sweep `{0.0, 0.1, ..., 0.9}`.
- Part 3: Q sweep `{5, 10, ..., 50}` (`p_ar = M*Q` varies, `p_white = p_total - p_ar`).
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}` (`p_ar`, `p_white`, and `s = M` all vary).

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_da/demo_trex_da_05_bt_ar1_block_white/demo_trex_da_05_bt_ar1_block_white.py
```

The file inserts both its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` and `da_sim_common` modules import directly. Output is
console-only; the folder carries its own `simulation_results/` directory for
optionally recorded summaries. Monte Carlo parts run in parallel via Python
`multiprocessing` using the `spawn` start method (up to `min(6, os.cpu_count())`
workers).

## Interpretation

- With most of the design being uncorrelated white noise, there is much less correlated false-positive-prone structure to latch onto, so FDR tends to be even more tightly controlled than Demo 04's all-block design at comparable SNR.
- Compare directly against Demo 04 (no white noise, same block AR(1) core) to isolate the effect of a much larger, mostly-null p versus a smaller fully-block design; the Q/M/linkage sweeps remain the key recovery diagnostics.

**Last updated**: 2026-07-08
