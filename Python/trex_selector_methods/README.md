# T-Rex Selector Methods — Python Demos

Python demos for the T-Rex selector family, using the `trex_selector_neo`
package. Each variant lives in its own sub-directory, following the C++
reference tree's **folder-per-demo** layout (one subfolder per demo, each with
its own `README.md` and `simulation_results/`, plus suite-level shared modules).

| Variant | Folder | Status |
|---|---|---|
| Classical T-Rex | [trex/](trex/README.md) | Available — demos 01–07 (single runs, MC simulations, L-loop strategies, memory-mapped pipelines, scalability benchmark) |
| Dependency-Aware T-Rex | [trex_da/](trex_da/README.md) | Available — demos 01–08 (AR(1), equicorrelated/BT, NN/MA, AR(1)-block, heavy-tailed blocks, prior-groups) |
| Grouped Variable Selection T-Rex | [trex_gvs/](trex_gvs/README.md) | Available — demos 01–08 (Hastie, scattered, mixed/negative-trap blocks, heavy-tailed, AR(1), ARMA, block benchmark; EN / EN+AUG / IEN) |
| Screen-TRex (ultra-high-dim screening) | [trex_screening/](trex_screening/README.md) | Available — demos 01–06 (Ordinary/Bootstrap-CI screening, memory-mapped, dependency-aware correlated designs, Biobank "Algorithm 1", solver backends) |
| T-Rex Sparse PCA | [trex_spca/](trex_spca/README.md) | Available — demo 01 (MC comparison of T-Rex SPCA, 2 solvers × 2 modes, vs. ordinary/oracle PCA on a sparse 3-factor model) |

---

**Last updated**: 2026-07-08
