# Demo 05: T-Rex+GVS on Heavy-Tailed (Student-t(3)) Blocks

## Purpose

Stress-test T-Rex+GVS on non-Gaussian data. This demo reuses the same four-block active/trap geometry as Demo 03 (Mixed Blocks), but replaces every Gaussian component of the data-generating process — latent factors, within-block perturbations, and response noise — with scaled Student-t(3) variables, testing robustness to heavy-tailed departures from the Gaussian assumption.

---

## Data Generation Parameters (`make_t3_equicorr_blocks_dgp`)

$$
X_{ij} = \sqrt{\rho}\, Z_{i,g(j)} + \sqrt{1-\rho}\, E_{ij}, \qquad Z_{\cdot,k},\, E_{ij} \sim t(3)/\sqrt{3},
$$

where $t(3)/\sqrt{3}$ denotes a Student-t distribution with 3 degrees of freedom, rescaled to unit variance. Response noise $\varepsilon_i$ is also drawn from $t(3)/\sqrt{3}$, scaled by $\sigma_y$.

- Same 4-block geometry as Demo 03: sizes $\{20, 50, 80, 65\}$; **blocks 0–2 active** ($\beta_j = 3$), **block 3 an inactive heavy-tailed trap**.
- $n = 200$, $p = 500$, $s = 150$.

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

- **EN** (`TENET`), **EN+AUG** (`TENET_AUG`), **IEN**.

---

## Three Parts

1. **SNR sweep** at fixed $\rho = 0.75$: $\mathrm{SNR} \in \{0.1, 0.2, 0.5, 1.0, 2.0, 5.0\}$.
2. **$\rho$ sweep** at fixed $\mathrm{SNR}=2.0$.
3. **2-D SNR $\times$ $\rho$ grid**.

---

## Output Files

Written to `simulation_results/`:

- `gvs_t3_blocks_snr.txt` / `.csv`
- `gvs_t3_blocks_rho.txt` / `.csv`
- `gvs_t3_blocks_2d.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_gvs/demo_trex_gvs_05_mc_sim_t3_blocks/demo_trex_gvs_05_mc_sim_t3_blocks
```

---

## Interpretation

- **EN** and **EN+AUG** are expected to remain reasonably well FDR-controlled with TPR rising toward the mid-$0.9$s as SNR grows, though slightly below the equivalent equicorrelated-Gaussian case (Demo 03) at comparable SNR — the heavy tails inject occasional large outliers that make exact group boundaries noisier to detect.
- Unlike the uniformly-positive-correlation Gaussian scenarios (Demos 01–03), **IEN** is expected to perform much less conservatively here, with FDR tightening (dropping well below the target as SNR grows) while TPR climbs into the $0.7$–$0.8$ range at higher SNR — a qualitatively different, more competitive pattern than IEN's behavior on the Gaussian equicorrelated blocks.
- The main question to track across the SNR sweep: does the heavy-tailed noise mainly cost power (lower TPR) or does it also compromise FDR control (values drifting above $0.1$)?

---

**Last updated**: 2026-07-08
