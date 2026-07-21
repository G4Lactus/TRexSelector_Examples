# Demo 06: Memory-Mapped Data — Single Run (Python)

## Purpose

Demonstrate memory-mapped (mmap) I/O workflows for the T-Rex selector via two
single-run patterns. Setting `use_memory_mapping = True` activates both
mmap storage of the internal dummy matrices `D` and solver LARS-path
serialization to disk between T-loop iterations; Part B additionally
memory-maps the design matrix `X` itself. Python port of
`cpp/trex_selector_methods/trex/demo_trex_06_mmap`.

## DGP / Data

Each part runs twice (low- then high-dimensional) via `dgp_gauss_snr`:

- Low-dimensional (`n > p`): `n = 5000`, `p = 1000`
- High-dimensional (`p > n`): `n = 1000`, `p = 5000`
- Fixed true support (0-based): `{27, 149, 398, 420, 4}` (`s = 5`)
- Fixed coefficients `beta_j = 1`, SNR = 1.0, tFDR = 0.1, seed = 58

Control (`_MMAP_CTRL_DICT`): TLARS, `K = 20`, `max_dummy_multiplier = 10`,
`lloop_strategy = "HCONCAT"`, `tloop_stagnation_stop = True`,
`use_memory_mapping = True`.

## Parts

- **Part A** (`part_a_d_mmap_solver_serial`): `X` stays in RAM;
  `use_memory_mapping = True` gives D-mmap + solver serialization.
- **Part B** (`part_b_full_mmap`): `X` is written to a temporary binary file via
  `numpy_to_memmap`, and the zero-copy `X_mmap.to_numpy()` view is passed to
  `TRexSelector`. The temp file is removed in a `finally` block (exception-safe,
  mirrors the C++ RAII `MmapFileGuard`).

`main()` runs Part A (low then high dim) followed by Part B (low then high dim).
Each run prints (and writes) the C++-style selection block via
`print_selection_results`: selected indices, FDP/TPP, sparse `phi_prime` /
`phi_mat`, `fdp_hat_mat`, `r_mat`, and the voting grid, plus `T_stop` and `L`.

## Imports

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `dgp_gauss_snr`, `trex_helpers`, and `trex_sim_common`
modules. The selector layer is used via the root alias
`import trex_selector_neo as tsn` (`tsn.TRexSelector`); `numpy_to_memmap`
comes from `trex_selector_neo.utils` via an explicit `from`-import.

## Output

Dual console + file output (mirrors the C++ `print_dual` pattern). Per
configuration, one `.txt` (full selection block) and one per-variable
selection `.csv` (via `save_selection_csv`) in this demo's own
`simulation_results/data/` folder, stems `d02_trex_mmap_demo_a_n{n}_p{p}` and
`d02_trex_mmap_demo_b_n{n}_p{p}` — matching the C++ result files. This is a
single-run demo (no `multiprocessing`).

## Running

```bash
python Python/trex_selector_methods/trex/demo_trex_06_mmap/demo_trex_06_mmap.py
```

**Last updated**: 2026-07-21
