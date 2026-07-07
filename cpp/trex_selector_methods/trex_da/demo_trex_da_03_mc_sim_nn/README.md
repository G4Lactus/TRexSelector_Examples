# Demo 03: DA-TRex on Banded (Nearest-Neighbor) Data

## Purpose

Study **`DAMethod::NN`**, which targets a banded (MA($\kappa$)) covariance structure where only $\kappa$ nearest-neighboring columns are correlated — as opposed to AR(1)'s geometrically-decaying correlation across *all* column distances.

> **Status**: All four parts run in `main()` and write output when executed. The committed checkout does not yet include the result files under `simulation_results/`.

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

## Output Files

Written to `simulation_results/` when run (7 scenario stems, one `.txt`+`.csv` pair each = 14 files):

- Part 1 (SNR): `da_trex_mc_da_nn_snr_capped.txt` / `.csv`, `da_trex_mc_da_nn_snr_random.txt` / `.csv`
- Part 2 ($\rho$): `da_trex_mc_da_nn_rho_capped.txt` / `.csv`, `da_trex_mc_da_nn_rho_random.txt` / `.csv`
- Part 3 ($\kappa$): `da_trex_mc_da_nn_kappa_capped.txt` / `.csv`, `da_trex_mc_da_nn_kappa_random.txt` / `.csv`
- Part 4 (2D $\kappa\times\rho$): `da_trex_mc_da_nn_kappa_rho.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_da/demo_trex_da_03_mc_sim_nn/demo_trex_da_03_mc_sim_nn
```

---

## Interpretation

- Since NN correlation has a hard cutoff at distance $\kappa$, expect DA-NN's correction to be most effective (best FDR control with least TPR cost) when active variables are spaced further apart than $\kappa$ — directly analogous to Demo 01's gap-vs-$\kappa$-boundary story, but with a sharp rather than gradually-decaying correlation profile.
- The $\kappa$ sweep (Part 3) and the 2D $\kappa\times\rho$ sweep (Part 4) are the most informative views here: as $\kappa$ grows, the effective correlation neighborhood widens, and DA-NN should need progressively wider correction windows (with a corresponding TPR cost analogous to Demo 01's findings for AR(1)).
- Compare this demo directly with **Demo 03b**, which applies the *same* NN correction to *AR(1)* data — Demo 03 is the "correctly specified" counterpart to that mismatch study.

---

**Last updated**: 2026-07-08
