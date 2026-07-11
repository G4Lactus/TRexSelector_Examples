# Demo 07: Monte Carlo Simulation with Memory-Mapped Data (Python)

## Purpose

Validate T-Rex selector FDR/TPR under memory-mapped (mmap) I/O over many Monte
Carlo trials, mirroring the single-run patterns of Demo 06. Python port of
`cpp/trex_selector_methods/trex/demo_trex_07_mc_sim_mmap`.

## DGP / Data

Both parts use the high-dimensional configuration:

- `n = 300`, `p = 1000`, `s = 10`
- True support drawn **once** via `make_support_random(s, p, seed=24)`, shared
  across all MC trials
- Fixed coefficients `beta_j = 1`, TLARS only
- SNR sweep: `{0.1, 0.5, 1.0, 2.0, 5.0}`
- `num_MC = 500` (`_NUM_MC = 500`)

Base control (`_MMAP_BASE_CTRL`): `K = 20`, `max_dummy_multiplier = 10`,
`lloop_strategy = "HCONCAT"`, `tloop_stagnation_stop = True`,
`tloop_max_stagnant_steps = 7`, `opt_threshold = 0.75`,
`use_memory_mapping = True`.

## Parts

- **Part C** (`trx_03_part_c`): in-memory `X` + `use_memory_mapping = True`
  (D-mmap + solver serialization), via `run_mc_trex`.
- **Part D** (`trx_03_part_d`): fully memory-mapped pipeline — each MC trial
  writes `X` to a unique temp binary file (`numpy_to_memmap`), passes
  `X_mmap.to_numpy()` to `TRexSelector`, and deletes the file in a `finally`
  block — via `run_mc_trex_full_mmap`. Adds X-mmap on top of D-mmap + solver
  serialization.

Both record mean FDR, TPR, Avg L, and Avg T per SNR (`track_L_T=True`).
`main()` runs Part C then Part D.

## Imports

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `support_generators` and `trex_sim_common` modules. The
package imports as `trex_selector_neo`.

## Output

Writes results to this demo's own `simulation_results/` folder via
`save_mc_results`, with stems `d03_trex_mmap_demo_c_n300_p1000_sw7` and
`d03_trex_mmap_demo_d_n300_p1000_sw7` (aligned `.txt` table plus tidy `.csv`).

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex/demo_trex_07_mc_sim_mmap/demo_trex_07_mc_sim_mmap.py
```

The MC loop is parallelized with Python `multiprocessing` (spawn start method)
using `_NUM_WORKERS = 6`.

**Last updated**: 2026-07-08
