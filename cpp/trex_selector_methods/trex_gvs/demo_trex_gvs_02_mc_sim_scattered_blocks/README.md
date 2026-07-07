# Demo 02: T-Rex+GVS on Scattered Blocks

## Purpose

Test whether spatially **scattering** the members of a correlated group across the design matrix — rather than keeping them contiguous — affects T-Rex+GVS group recovery. This isolates the effect of column layout from the effect of correlation strength, since the underlying group structure is otherwise identical to Demo 01.

---

## Data Generation Parameters (`make_scattered_grouped_dgp`)

$$
X_{ij} = Z_{i,\,g(j)} + \sigma_x\, \xi_{ij}, \qquad \xi_{ij} \sim \mathcal{N}(0,1),
$$

- `n_groups` correlated groups of `group_size` variables each; here 3 groups of 50 variables, matching Demo 01's Hastie design.
- Group membership is assigned via a **random column permutation** instead of contiguous blocks — group members can be spread arbitrarily across the $p$ columns.
- All grouped variables are active, $\beta_j = 3$; background columns are i.i.d. $\mathcal{N}(0,1)$, $\beta_j = 0$.
- $n = 200$, $p = 500$, $s = 150$.
- Within-group correlation $\rho = 1/(1+\sigma_x^2)$, same convention as Demo 01.

---

## Control Parameters

```
K = 20
tFDR = 0.1
corr_max = 0.98
hc_linkage = Single
lambda2_method = CV_1SE_CCD
MC = 200
```

---

## Methods Compared

- **EN** (`TENET`), **EN+AUG** (`TENET_AUG`), **IEN** — same three methods as Demo 01.

---

## Three Parts

1. **SNR sweep** at fixed $\sigma_x = 0.1$: $\mathrm{SNR} \in \{0.1, 0.2, 0.5, 1.0, 2.0, 5.0\}$.
2. **$\rho$ sweep** at fixed $\mathrm{SNR} = 2.0$.
3. **2-D SNR $\times$ $\rho$ grid**.

---

## Output Files

Written to `simulation_results/`:

- `gvs_scattered_grouped_snr.txt` / `.csv`
- `gvs_scattered_grouped_rho.txt` / `.csv`
- `gvs_scattered_grouped_2d.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_gvs/demo_trex_gvs_02_mc_sim_scattered_blocks/demo_trex_gvs_02_mc_sim_scattered_blocks
```

---

## Interpretation

- Because HAC-based grouping clusters columns by **correlation**, not by column adjacency, scattering group members across the design is expected to leave FDR/TPR patterns very close to Demo 01's contiguous-block results — the clustering step is insensitive to column order by construction.
- This demo is best read side-by-side with Demo 01: similar results confirm that group *discovery* in this pipeline depends on correlation structure rather than spatial layout, while a noticeable divergence would point to a layout-sensitive bug or clustering limitation.
- As in Demo 01, expect **IEN** to be conservative (very low FDR and TPR) on this all-active-groups design.

---

**Last updated**: 2026-07-08
