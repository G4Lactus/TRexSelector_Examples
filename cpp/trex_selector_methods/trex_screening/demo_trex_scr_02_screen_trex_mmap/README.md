# Demo 02: Screen-TRex (Memory-Mapped)

## Purpose

Repeat Demo 01's Ordinary vs. Bootstrap-CI comparison with **memory-mapped** dummy-matrix storage (`use_memory_mapping = true`), to validate that Screen-TRex's disk-backed workflow produces the same behavior as the in-memory version while enabling larger, RAM-constrained problems.

> **Status**: `simulation_results/` is currently empty — this demo has not yet been run in this checkout.

---

## Data Generation Parameters

Identical to Demo 01:

- Single-run demos: $n=300$, $p=1000$, fixed support $\{27, 149, 398, 4, 42\}$ ($s=5$), $\beta_j=1.0$, $\mathrm{SNR}=5.0$.
- MC sweep: $n=300$, $p=1000$, random $s=10$ support per trial, $\mathrm{SNR} \in \{0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0\}$, $\mathrm{MC}=500$.

---

## Control Parameters

```
K = 20
solver = TLARS
use_memory_mapping = true     # only difference from Demo 01
R_boot = 1000
ci_grid_step = 0.001
```

---

## Methods Compared

- **Screen-TRex Ordinary (mmap)**
- **Screen-TRex Bootstrap-CI (mmap)**

---

## Output Files (expected)

Written to `simulation_results/` once run:

- `d02_screen_trex_mmap_mc_n300_p1000.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_scr_02_screen_trex_mmap
```

---

## Interpretation

- This demo is best read side-by-side with Demo 01: since the DGP, control parameters, and methods are identical, the FDR/TPR/estimated-FDR curves are expected to match Demo 01 closely (up to Monte Carlo noise from different random seeds) — the only intended difference is where the dummy matrices are stored (disk vs. RAM), not the statistical behavior of the selector.
- A material divergence from Demo 01's results at the same SNR points would suggest an mmap-specific issue (e.g. a stale or improperly synchronized backing file) rather than a genuine statistical effect.

---

**Last updated**: 2026-07-04
