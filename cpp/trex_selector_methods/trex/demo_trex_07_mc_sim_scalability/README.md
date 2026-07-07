# Demo 07: Scalability Benchmark (Placeholder)

## Status

**Work in progress**: This demo folder is reserved for a future large-scale scalability benchmark. The source file `demo_trex_07_mc_sim_scalability.cpp` is currently an **empty placeholder**, and its `add_demo_target(...)` entry is **commented out** in `../CMakeLists.txt`, so nothing is built or run for this demo yet. Everything below is a proposed design, not implemented behaviour.

---

## Proposed Purpose

Characterize **T-Rex selector runtime and memory usage** as $n$ and $p$ scale to extreme sizes. This would include:

- **Runtime scaling**: Measure wall-clock time as a function of $n$, $p$, and $K$ (number of random experiments)
- **Memory scaling**: Track peak RAM usage and mmap file sizes
- **Solver comparison**: Compare computational efficiency across different LARS variants
- **Parallel scaling**: Quantify speedup from OpenMP parallelization (if enabled)

---

## Potential Design

### Scenarios
| Scenario | $n$ | $p$ | Storage | Purpose |
|----------|-----|-----|---------|---------|
| A | 10,000 | 10,000 | In-memory | Moderate scale, RAM-resident |
| B | 10,000 | 100,000 | Mmap | High-dimensional, disk I/O |
| C | 100,000 | 10,000 | Mmap | Large sample, high-dimensional mix |
| D | 100,000 | 100,000 | Mmap | Extreme scale (if hardware permits) |

### Output
- Timing breakdown (data generation, L-loop, T-loop per iteration, total)
- Memory peak (RSS, mmap file size)
- Scaling plots (runtime vs. $p$, memory vs. $n \times p$)
- Solver performance comparison

---

## To Do

- [ ] Design and implement large-scale DGP
- [ ] Set up timing instrumentation for each algorithm phase
- [ ] Implement memory profiling hooks
- [ ] Configure appropriate $K$ for different scales (trade-off: larger $K$ → better FDP estimates, higher cost)
- [ ] Run and analyze scalability results
- [ ] Create visualization and summary plots

---

## Placeholder

When implemented, this demo will demonstrate the **computational efficiency and scalability** of the T-Rex selector for production-scale variable selection tasks.

---

**Last updated**: 2026-07-08  
**Status**: Planned, not yet implemented (empty source, excluded from the build)
