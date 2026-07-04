# Demo 03: DA-TRex on Banded (Nearest-Neighbor) Data

## Purpose

Study **`DAMethod::NN`**, which targets a banded (MA($\kappa$)) covariance structure where only $\kappa$ nearest-neighboring columns are correlated — as opposed to AR(1)'s geometrically-decaying correlation across *all* column distances.

> **Status**: `simulation_results/` is currently empty — this demo has not yet been run in this checkout.

---

## Data Generation Parameters (`dgp_nn`)

$$
X_{i,j} = \sum_{l=0}^{\kappa} \theta_l\, \eta_{i,j+l}, \qquad \theta_l \propto \rho^l \text{ (normalized to unit variance)}.
$$

- $n=300$, $p=1000$, `kappa=3` (bandwidth), `rho=0.7` (MA coefficient decay).
- Unlike AR(1), correlation between columns $j$ and $k$ is **exactly zero** once $|j-k|>\kappa$ — a genuinely banded, finite-range dependency.

---

## Four Parts

1. **SNR sweep** — $\kappa=3$, $\rho=0.7$ fixed; `CappedSpread` and `Random` support.
2. **$\rho$ sweep** — $\kappa=3$, $\mathrm{SNR}=2.0$ fixed; dual support.
3. **$\kappa$ sweep** — $\rho=0.7$, $\mathrm{SNR}=2.0$ fixed; dual support.
4. **2D $\kappa \times \rho$ sweep** — $\mathrm{SNR}=2.0$ fixed; all 3 solvers; `CappedSpread` and `Random` sub-sections.

All parts: $n=300$, $p=1000$, $s=10$, $K=20$, $\mathrm{tFDR}=0.2$, $\mathrm{MC}=200$, `base_seed=2026`.

---

## Output Files (expected)

Written to `simulation_results/` once run, following the shared `da_trex_mc_{scenario_tag}.{txt,csv}` naming convention (verify exact tag strings in the source if you rerun this demo).

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_da_03_mc_sim_nn
```

---

## Interpretation

- Since NN correlation has a hard cutoff at distance $\kappa$, expect DA-NN's correction to be most effective (best FDR control with least TPR cost) when active variables are spaced further apart than $\kappa$ — directly analogous to Demo 01's gap-vs-$\kappa$-boundary story, but with a sharp rather than gradually-decaying correlation profile.
- The $\kappa$ sweep (Part 3) and the 2D $\kappa\times\rho$ sweep (Part 4) are the most informative views here: as $\kappa$ grows, the effective correlation neighborhood widens, and DA-NN should need progressively wider correction windows (with a corresponding TPR cost analogous to Demo 01's findings for AR(1)).
- Compare this demo directly with **Demo 03b**, which applies the *same* NN correction to *AR(1)* data — Demo 03 is the "correctly specified" counterpart to that mismatch study.

---

**Last updated**: 2026-07-04
