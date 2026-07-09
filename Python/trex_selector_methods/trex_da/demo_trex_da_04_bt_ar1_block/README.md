# Demo 04: T-Rex+DA+BT on Block-Diagonal AR(1) Data (Python)

## Purpose

Introduce the `BT` DA method (binary-tree / dendrogram HAC clustering) on a clean
block-diagonal AR(1) design, and study sensitivity to the HAC linkage method
(single / complete / average) across four sweep dimensions: SNR, rho, block size Q,
and number of blocks M.

Python port of `R/trex_selector_methods/trex_da/demo_trex_da_04_bt_ar1_block` and
`cpp/.../demo_trex_da_04_mc_sim_bt_ar1_block`.

## Data-generating process

`dgp_ar1_block_snr` partitions p predictors into M independent blocks of size Q
(`p = M * Q`), with an AR(1) recursion within each block and zero correlation between
blocks. One active variable per block, so `s = M` (`make_support_bt_one_per_block`).
Base config: `n=150, M=5, Q=5` (p=25, s=5), `rho=0.7, snr=2.0, amplitude=1.0,
tFDR=0.2, K=20, num_MC=200, seed=2026`.

## DA method

`build_da_selector("BT", ..., hc_linkage=lnk)` constructs a
`trex_selector_neo.TRexDASelector` with a `TRexDAControlParameter`
(`method = DAMethod.BT`, `hc_linkage` set from `LinkageMethod.{Single,Complete,Average}`)
and a `TRexControlParameter` (`K = 20`, `solver_type = SolverTypeForTRex.TLARS`). Every
sweep loops the linkage over `{single, complete, average}`.

## What it runs

All parts active by default (`RUN_PART_1..4 = True`). Each is heavy
(num_MC x sweep points x 3 linkages) — reduce `num_MC` to preview.

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- Part 2: rho sweep `{0.0, 0.1, ..., 0.9}`.
- Part 3: Q sweep `{5, 10, ..., 50}` (block size; `p = M*Q` varies, s=M=5 fixed).
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}` (number of blocks; `p = M*Q` and `s = M` both vary).

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_da/demo_trex_da_04_bt_ar1_block/demo_trex_da_04_bt_ar1_block.py
```

The file inserts both its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` and `da_sim_common` modules import directly. Output is
console-only; the folder carries its own `simulation_results/` directory for
optionally recorded summaries. Monte Carlo parts run in parallel via Python
`multiprocessing` using the `spawn` start method (up to `min(6, os.cpu_count())`
workers).

## Interpretation

- FDR/TPR patterns depend on how well the HAC dendrogram cut recovers the true block structure; recovering every active block consistently is harder than recovering isolated variables, so TPR tends to plateau below 1.0 even at high SNR.
- The Q and M sweeps are the key diagnostics for how block size and block count affect recovery; the linkage sweep answers whether single/complete/average materially changes FDR/TPR on this structure.

**Last updated**: 2026-07-08
