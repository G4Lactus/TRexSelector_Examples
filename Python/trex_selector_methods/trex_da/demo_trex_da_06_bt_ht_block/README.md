# Demo 06: T-Rex+DA+BT on Heavy-Tailed Block Data (Python)

## Purpose

Stress-test `BT`-based DA-TRex under heavy-tailed (Student-t) predictors and/or
noise on a block-diagonal design. Studies both statistical robustness (Gaussian vs
heavy-tailed scenarios) and sensitivity across five sweep dimensions, plus a
dedicated linkage-method sweep.

Python port of `R/trex_selector_methods/trex_da/demo_trex_da_06_bt_ht_block` and
`cpp/.../demo_trex_da_06_mc_sim_bt_ht_block_sweeps`. (The C++ variant uses a
Cholesky heavy-tail DGP; the intent matches.)

## Data-generating process

`dgp_ht_block_snr` draws AR(1)-block Gaussian predictors and converts them to a
multivariate Student-t(nu) via a scale mixture; one active variable per block
(`s = M`). Two noise scenarios are run for every sweep: Scenario 1 uses Gaussian
noise (`heavy_noise=False`), Scenario 2 uses Student-t(nu) noise
(`heavy_noise=True`). Base config: `n=150, M=5, Q=5` (p=25, s=5), `rho=0.8, nu=3.0,
snr=2.0, amplitude=1.0, tFDR=0.2, K=20, num_MC=200, seed=2026`.

## DA method

`build_da_selector("BT", ..., hc_linkage=lnk)` constructs a
`trex_selector_neo.TRexDASelector` with a `TRexDAControlParameter`
(`method = DAMethod.BT`, `hc_linkage` from `LinkageMethod.{Single,Complete,Average}`)
and a `TRexControlParameter` (`K = 20`, `solver_type = SolverTypeForTRex.TLARS`). Each
non-linkage sweep loops the linkage over `{single, complete, average}` inside each
noise scenario.

## What it runs

All parts active by default (`RUN_PART_1..6 = True`). Each is heavy
(2 scenarios x 3 linkages x sweep points x num_MC) — reduce `num_MC` to preview.

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`.
- Part 2: rho sweep `{0.0, 0.1, ..., 0.9}`.
- Part 3: Q sweep `{5, 10, ..., 50}`.
- Part 4: M sweep `{1, 3, 5, 10, 15, 20, 30}`.
- Part 5: tFDR sweep `{0.05, 0.10, ..., 0.50}`.
- Part 6: linkage sweep (linkage is the swept parameter), per scenario.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_da/demo_trex_da_06_bt_ht_block/demo_trex_da_06_bt_ht_block.py
```

The file inserts both its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` and `da_sim_common` modules import directly. Output is
console-only; the folder carries its own `simulation_results/` directory for
optionally recorded summaries. Monte Carlo parts run in parallel via Python
`multiprocessing` using the `spawn` start method (up to `min(6, os.cpu_count())`
workers).

## Interpretation

- Under heavy-tailed noise, FDR is typically less tightly controlled (hovering around or above the tFDR=0.2 target) and TPR plateaus lower than the Gaussian block demos, reflecting the added estimation difficulty from occasional large outliers in t(3) data.
- Compare the Gaussian vs heavy scenarios (same base parameters, only the noise/predictor distribution differs) to isolate the cost of heavy tails; with noisier correlation estimates, linkage-method choice may matter more here than in the cleaner Gaussian demos. The tFDR sweep checks whether realized FDR tracks the target across levels.

**Last updated**: 2026-07-08
