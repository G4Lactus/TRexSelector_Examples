# Demo 07: T-Rex+DA+BT on Heavy-Tailed Block + White-Noise Data (Python)

## Purpose

Combine Demo 05's white-noise dilution with Demo 06's heavy-tailed robustness study:
heavy-tailed AR(1)-Toeplitz blocks embedded in a larger design padded with i.i.d.
Student-t white-noise columns. Sensitivity across five sweep dimensions under two
noise scenarios.

Python port of `R/trex_selector_methods/trex_da/demo_trex_da_07_bt_ht_block_white`
and `cpp/.../demo_trex_da_07_mc_sim_bt_ht_block_white`.

## Data-generating process

`dgp_ht_block_white_snr` builds M heavy-tailed Student-t(nu) AR(1) blocks of size Q
(`p_ar = M * Q`, the actives) plus `p_white = p_total - p_ar` heavy-tailed white-noise
columns, with `p_total` fixed at 500. Active variables (`s = M`, one per AR block) lie
in the heavy-tailed AR part only. Two noise scenarios per sweep: Scenario 1 Gaussian
noise (`heavy_noise=False`), Scenario 2 Student-t(nu) noise (`heavy_noise=True`). Base
config: `n=150, p_total=500, M=5, Q=5` (p_ar=25, p_white=475, s=5), `rho=0.8, nu=3.0,
snr=2.0, amplitude=1.0, tFDR=0.2, K=20, num_MC=200, seed=2026`.

## DA method

`build_da_selector("BT", ..., hc_linkage=lnk)` constructs a
`trex_selector_neo.TRexDASelector` with a `TRexDAControlParameter`
(`method = DAMethod.BT`, `hc_linkage` from `LinkageMethod.{Single,Complete,Average}`)
and a `TRexControlParameter` (`K = 20`, `solver_type = SolverTypeForTRex.TLARS`). Every
sweep loops the linkage over `{single, complete, average}` inside each noise scenario.

## What it runs

All parts active by default (`RUN_PART_1..5 = True`). There is no dedicated
linkage-only section (unlike Demo 06) — linkage is swept as an outer loop within each
part. Each part is heavy — reduce `num_MC` to preview.

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- Part 2: rho sweep `{0.0, 0.1, ..., 0.9}`.
- Part 3: Q sweep `{5, 10, ..., 50}` (`p_ar = M*Q` varies, `p_white = p_total - p_ar`).
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}` (`p_ar`, `p_white`, and `s = M` all vary).
- Part 5: tFDR sweep `{0.05, 0.10, ..., 0.50}`.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_da/demo_trex_da_07_bt_ht_block_white/demo_trex_da_07_bt_ht_block_white.py
```

The file inserts both its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` and `da_sim_common` modules import directly. Output is
console-only; the folder carries its own `simulation_results/` directory for
optionally recorded summaries. Monte Carlo parts run in parallel via Python
`multiprocessing` using the `spawn` start method (up to `min(6, os.cpu_count())`
workers).

## Interpretation

- Read this as the "combined effects" check: compare against Demo 04 (Gaussian, no white noise), Demo 05 (Gaussian, with white noise), and Demo 06 (heavy-tailed, no white noise) to decompose how heavy tails versus white-noise dilution each affect FDR/TPR, and whether their effects compound or partially offset.
- As with the other block demos, the Q/M/linkage sweeps probe recovery under the white-noise-diluted design, and the tFDR sweep checks whether realized FDR tracks the target across levels.

**Last updated**: 2026-07-08
