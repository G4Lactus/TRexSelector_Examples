# Demo 08: T-Rex+GVS Block Benchmark — HAC-Discovered vs. Oracle Groups

## Purpose

Benchmark how much accurate **group discovery** (via HAC clustering) costs compared to knowing the true groups in advance ("oracle" groups), across three different within-block truth patterns. This is the most structurally different demo in the folder: instead of an SNR/$\rho$ sweep for a single method set, it runs a **4-method $\times$ 3-scenario** comparison over a $(\rho, \mathrm{SNR})$ grid.

> **Note**: this demo prints all results to the console only — its `simulation_results/` folder is intentionally empty (no `.txt`/`.csv` files are written).

---

## Data Generation Parameters (`make_block_equicorr_dgp`)

$$
X_{ij} = \sqrt{\rho}\, Z_{i,g(j)} + \sqrt{1-\rho}\, \epsilon_{ij}, \qquad Z_{\cdot,k}, \epsilon_{ij} \sim \mathcal{N}(0,1).
$$

- $G = 100$ equal-size blocks of `block_size = 5` variables each ($p = 500$).
- `n_active_blocks = 10` blocks are chosen at random as active; each receives a random sign coefficient $b = \pm 3.0$.
- $n = 200$.

### Three truth scenarios (`BlockScenario`)

| Scenario | Active pattern within each active block | Resulting $s$ |
|---|---|---|
| **INDIVIDUAL** | Exactly 1 variable active per active block | $s = 10$ |
| **REPRESENTATIVE** | 2–3 random variables active per active block | $s \approx 20$–$30$ |
| **WHOLE** | All 5 variables in the block are active | $s = 50$ |

---

## Four Methods Compared

| Method | `GVSType` | Groups |
|---|---|---|
| **M1 (EN-C)** | EN | HAC-clustered (discovered from data via `corr_max` threshold, single-linkage) |
| **M2 (EN-O)** | EN | Oracle (true `prior_groups` from the DGP) |
| **M3 (IE-C)** | IEN | HAC-clustered |
| **M4 (IE-O)** | IEN | Oracle |

---

## Control Parameters

```
n = 200, p = 500, G = 100, block_size = 5, n_active_blocks = 10
tFDR = 0.1
K = 20
corr_max = 0.5
b = 3.0                      # nonzero coefficient magnitude
rho_grid = {0.5, 0.9}
snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}
MC = 500 per (rho, snr) cell, per scenario
```

Seeds are staggered per grid cell (`base_seed + cell_index * num_MC`) so that every $(\rho, \mathrm{SNR})$ cell across all three scenarios draws from a distinct, non-overlapping seed range.

---

## Metrics Reported

Beyond the usual coordinate-level FDR/TPR, this demo reports **block-level diagnostics**:

- **`blk_hit`** (block hit rate): share of active blocks with at least one selected variable (INDIVIDUAL/REPRESENTATIVE) — or **`full_blk`** (full-block rate): share of active blocks recovered in their entirety (WHOLE).
- **`block_FDR`**: share of hit blocks that are actually null.
- **`purity`**: block-purity rate — how consistently a true block maps to a single discovered group (always $1.0$ for oracle groups by construction; $<1.0$ for HAC-clustered groups indicates clustering did not perfectly separate blocks).
- **`T_stop`**, **`M_found`**: selector diagnostics (T-loop stopping iteration, number of groups found/selected).

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_gvs_08_mc_sim_block_bench
```

All three scenarios run in the same invocation and print one grid table per scenario, per $\rho$ value.

---

## Interpretation

- The central comparison is **M1 vs. M2** (and **M3 vs. M4**): if `corr_max` is well matched to the true within-block $\rho$, HAC-clustered groups should track oracle-group performance closely (`purity` near $1.0$, similar TPR/FDR); a persistent gap indicates the clustering threshold is splitting or merging true blocks incorrectly.
- The **INDIVIDUAL** scenario is expected to be the hardest (only 1 of 5 correlated variables per block is truly active, so a discovered group can "win" on a correlated null neighbor); **WHOLE** should be the easiest, since the entire block carries signal.
- Comparing $\rho = 0.5$ vs. $\rho = 0.9$ across scenarios shows how much within-block correlation strength helps or hurts group discovery — stronger correlation generally makes HAC clustering more accurate but can also make it harder to localize signal within a block in the INDIVIDUAL/REPRESENTATIVE scenarios.
- **EN vs. IEN** (M1/M2 vs. M3/M4) can be compared directly at matched group-source (clustered or oracle) to see whether the iterative elastic net changes the HAC-vs-oracle gap itself, or only the overall detection level.

---

**Last updated**: 2026-07-04
