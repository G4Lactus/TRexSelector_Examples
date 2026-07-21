# Demo 08: Scalability Benchmark (Python)

## Purpose

Characterize T-Rex selector **runtime and memory usage** as `n` and `p` scale,
across three base solvers (TLARS, TOMP, TAFS), on the **fully memory-mapped
pipeline**: every Monte Carlo trial generates `X` directly into an mmap
backing file (chunked column generation over `MemoryMappedMatrix`, so `X` is
never materialized in RAM), and `use_memory_mapping = True` additionally
memory-maps the internal dummy matrices `D` and serializes solver warm-start
state to disk. Python port of
`cpp/trex_selector_methods/trex/demo_trex_08_mc_sim_scalability`.

The benchmark reports, per scenario x solver x SNR grid point:

- **Runtime scaling**: wall-clock per phase (DGP, `select()`,
  per-`T`-iteration, total)
- **Memory scaling**: post-`select()` RSS, physical footprint, process peak
  RSS, and the `X`/`y` data sizes
- **Selection quality at scale**: FDR, TPR, Avg L, Avg T as a sanity check

It is a **2D sweep**: every scenario x solver combination runs over the SNR
grid `{0.5, 1.0, 2.0}`.

## Scenarios

| Scenario | n | p | num_MC | Purpose |
|---|---|---|---|---|
| A | 300 | 1,000 | 10 | Baseline |
| B | 1,000 | 10,000 | 10 | Moderate scale |
| C | 5,000 | 100,000 | 10 | Large sample x high-dimensional |

A tiny smoke-test scenario (`SCENARIOS_SMOKE`, `n = 150`, `p = 500`,
`num_MC = 2`) is gated behind `if False:` in `main()` for end-to-end pipeline
verification — mirrors the C++ `make_smoke_scenarios()` gate.

Scenario C is heavy: expect minutes per trial and tens of GiB of
swap/compression pressure from the per-experiment coefficient-path heap (see
the C++ README's *Memory behaviour* discussion — it applies unchanged, since
the selector core is the same C++ code).

## DGP / Data

Mirrors **Demo 03** (variable support):

- True support: `s = 10`, indices redrawn per trial via
  `make_support_random(seed=trial_seed)` (internal `+500_000` offset mirrors
  the C++ `draw_support`)
- Fixed coefficients `beta_j = 1` (`rnd_coef=False`)
- Seeding: `base_seed = 24 + 1000 * (3 * scenario_index + snr_index)`, unique
  per (scenario, SNR) and shared across solvers so every solver sees identical
  data; the selector runs with `seed=-1` (hardware entropy)
- Trials run **sequentially** (no `ProcessPoolExecutor`) so wall-clock and RSS
  readings are undistorted; the random experiments inside `select()` may still
  run in parallel (`parallel_rnd_experiments` keeps its default, mirroring the
  C++ OpenMP setup)

## Control Parameters

```text
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Max dummies L = 10p
use_max_T_stop = True            # Cap T <= ceil(n/2)
lloop_strategy = STANDARD        # Fresh i.i.d. dummy matrix per L-loop iteration
tloop_stagnation_stop = True     # Early exit when R_mat stagnates
tloop_max_stagnant_steps = 5     # Stagnation window
use_memory_mapping = True        # All scenarios: D mmap + solver serialization
tFDR = 0.1                       # Target FDR control level
```

## Measurements

Twelve metric rows per solver (mirrors the C++ table): `FDR`, `TPR`, `Avg L`,
`Avg T`, `DGP s`, `Select s`, `s/T-iter`, `Total s`, `RSS MiB`, `FootprMiB`,
`PeakRSS`, `X+y MiB`.

Python adaptations of the memory hooks (cpp reads mach/`/proc` directly):
RSS via `ps -o rss=` (macOS) / `/proc/self/statm` (Linux); physical footprint
via `task_info(TASK_VM_INFO).phys_footprint` through `ctypes` (macOS) /
`VmRSS + VmSwap` (Linux), falling back to RSS on failure; peak RSS via
`resource.getrusage` (monotone over the process lifetime). `y` is kept in RAM
(`n` doubles — negligible); `X+y MiB` counts the `X` backing file plus `y`'s
in-memory bytes.

## Imports

The demo inserts its own folder and the parent suite dir onto `sys.path` to
import the shared `support_generators`, `trex_helpers`, and `trex_sim_common`
modules. The selector layer is used via the root alias
`import trex_selector_neo as tsn` (`tsn.TRexSelector`); `MemoryMappedMatrix`
and `AccessMode` come from `trex_selector_neo.utils` via explicit
`from`-imports.

## Output Files

Written to `simulation_results/data/` with stem
`demo_trex_08_scalability_results`, **re-saved after each completed scenario**
(a crash or OOM kill on a later scenario preserves the finished ones):

- `demo_trex_08_scalability_results.txt` — scenario info block plus an aligned
  table in the demo 02–07 style, one column per scenario x SNR combination
  (labelled `<scenario>/<snr>`, e.g. `A/0.5`), twelve metric rows per solver.
- `demo_trex_08_scalability_results.csv` — tidy long format, header column
  order `scenario,solver,n,p,num_mc,snr,metric,value` (same schema the C++
  suite plotter's `scalability` mode reads).

## Running

```bash
python Python/trex_selector_methods/trex/demo_trex_08_mc_sim_scalability/demo_trex_08_mc_sim_scalability.py
```

Scenarios run outermost (every solver finishes the current scenario before the
next, larger one starts). The full A/B/C benchmark is long — scenario C
dominates; flip the `if False:` smoke gate in `main()` for a quick end-to-end
check.

**Last updated**: 2026-07-21
