# Demo 05: Biobank Screen-TRex (Memory-Mapped)

## Purpose

Repeat Demo 04's Biobank Screen-TRex ("Algorithm 1") multi-phenotype workflow with **memory-mapped** dummy-matrix storage (`use_memory_mapping = true`), validating that the adaptive Ordinary → Bootstrap-CI → T-Rex-fallback routing behaves the same way when D-matrices are disk-backed — relevant for biobank-scale problems where in-memory storage becomes infeasible.

> **Status**: `simulation_results/` is currently empty — this demo has not yet been run in this checkout.

---

## Data Generation Parameters

Identical to Demo 04:

- Single phenotype: $n=300$, $p=1000$, support $\{4,27,42,149,398\}$ ($s=5$), $\mathrm{SNR}=1.0$.
- Multiple phenotypes: $q=5$, shared $\boldsymbol{X}$, per-phenotype support sizes $\{5,6,5,7,5\}$, $\mathrm{SNR}\in\{1,2,5,10,20\}$.
- MC single phenotype: $s=10$, $\mathrm{SNR}\in\{0.1,0.5,1.0,2.0,5.0\}$, $\mathrm{MC}=500$.
- MC multiple phenotypes: $q=5$, $s=5$ per phenotype, same SNR grid, $\mathrm{MC}=500$.

---

## Control Parameters

```
lower_bound_FDR = 0.05
upper_bound_FDR = 0.15
target_FDR_trex = 0.10
R_boot = 1000
ci_grid_step = 0.001
K = 20
solver = TLARS
use_memory_mapping = true     # only difference from Demo 04
```

---

## Methods (Adaptive Routing)

Same three-tier routing as Demo 04: **Screen-TRex (Ordinary)** → **Screen-TRex (Bootstrap-CI)** → **T-Rex (fallback)**.

---

## Output Files (expected)

Written to `simulation_results/` once run:

- `d05_biobank_screen_trex_mc_snr_n300_p1000_s10_mmap.txt` / `.csv` — single-phenotype MC SNR sweep.
- `d05_biobank_screen_trex_mc_multi_n300_p1000_q5_s5_mmap.txt` / `.csv` — multi-phenotype MC SNR sweep.

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_scr_05_biobank_screen_trex_mmap
```

---

## Interpretation

- Read alongside Demo 04: since the DGP, routing logic, and control parameters are identical, Usage %, conditional FDP/TPP, and estimated FDR are all expected to closely match Demo 04's results at the same SNR points — the only intended difference is D-matrix storage location.
- This demo is primarily a feasibility/parity check: confirming that the adaptive routing logic does not depend on assumptions that break under memory-mapped storage (e.g. thread-local temporary file handling across the many phenotype × trial combinations).

---

**Last updated**: 2026-07-04
