# Demo 05: Memory-Mapped Data — Single Run

## Purpose

Demonstrate memory-mapped (mmap) I/O workflows for the T-Rex selector via two
single-run usage patterns. `use_memory_mapping = TRUE` memory-maps the internal dummy
matrices `D` and serializes solver LARS-path state to disk between T-loop iterations;
Part B additionally memory-maps the design matrix `X` itself.

## Data-generating process

`dgp_gauss_snr` — i.i.d. Gaussian design, `y = X beta + eps` calibrated to SNR.
Shared config (`DEMO_CFG`): fixed true support `c(28, 150, 399, 421, 5)` (1-based;
C++ 0-based `{27, 149, 398, 420, 4}`, `s = 5`), coefficients `beta_j = 1`, `SNR = 1.0`,
`tFDR = 0.1`, `seed = 58`. Each part runs twice: low-dimensional `n = 5000, p = 1000`
then high-dimensional `n = 1000, p = 5000`. Control via `make_mmap_trex_ctrl(solver = "TLARS")`.

## Parts

- Part A (`run_part_a`) — in-memory `X` + `use_memory_mapping = TRUE`: activates D-mmap
  and solver serialization inside `TRexSelector`. Mirror of the C++ Demo A.
- Part B (`run_part_b`) — memory-mapped `X` + `use_memory_mapping = TRUE`: `X` is written
  to a temporary mmap file via `execute_with_temp_mmap` before being passed to the
  selector (full pipeline: X mmap + D mmap + solver serialization). The R interface
  requires writing the generated matrix to disk first (C++ uses `SyntheticDataMapped`,
  where `X` is never fully in RAM).

`.print_and_save_d02` prints selected indices, FDP, TPP, `T_stop`, and `L`.

## Running

```bash
Rscript R/trex_selector_methods/trex/demo_trex_05_mmap/demo_trex_05_mmap.R
```

Each run writes a per-variable selection **CSV** (no `.txt`) to this demo's own
`simulation_results/` folder — header `variable_index,phi_prime,selected,true_positive`,
stems `d02_trex_mmap_demo_a_n{n}_p{p}.csv` and `d02_trex_mmap_demo_b_n{n}_p{p}.csv`. The
script sources `../../simulation_utils.R` and `../trex_sim_utils.R`.

**Last updated**: 2026-07-08
