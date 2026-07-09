# Demo 06: T-Rex+GVS on AR(1) Blocks

## Purpose

Evaluate T-Rex+GVS when within-block correlation follows an **AR(1) time-series structure** instead of equicorrelation. This uses the same four-block active/trap geometry as Demo 03, but replaces the shared-latent-factor equicorrelated model with an AR(1) Toeplitz covariance within each block, testing robustness to more realistic, decaying correlation patterns.

---

## Data Generation Parameters (`make_ar1_blocks_dgp`)

Within-block covariance follows an AR(1) Toeplitz structure,

$$
\Sigma_k(i,j) = \rho^{\lvert i-j \rvert},
$$

generated as $X_{\text{block}} = Z \cdot L^\top$, where $Z \sim \mathcal{N}(0, I_{n \times b_k})$ and $L$ is the Cholesky factor of $\Sigma_k$.

- Same 4-block geometry as Demo 03: sizes $\{20, 50, 80, 65\}$; **blocks 0–2 active** ($\beta_j = 3$), **block 3 an inactive AR(1) trap**.
- Gap columns (between/around blocks) are i.i.d. $\mathcal{N}(0,1)$.
- $n = 200$, $p = 500$, $s = 150$.

Unlike the equicorrelated model, correlation between two variables in the same AR(1) block decays with their column distance, so a block no longer behaves like a single tight cluster — variables at opposite ends of a large block can be nearly uncorrelated.

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

1. **SNR sweep** at fixed $\rho = 0.85$: $\mathrm{SNR} \in \{0.1, 0.2, 0.5, 1.0, 2.0, 5.0\}$.
2. **$\rho$ sweep** at fixed $\mathrm{SNR}=2.0$ (expanded coverage, 11 levels including $0.40$–$0.99$).
3. **2-D SNR $\times$ $\rho$ grid**.

---

## Output Files

Written to `simulation_results/`:

- `gvs_ar1_blocks_snr.txt` / `.csv`
- `gvs_ar1_blocks_rho.txt` / `.csv`
- `gvs_ar1_blocks_2d.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_gvs/demo_trex_gvs_06_mc_sim_ar1_blocks/demo_trex_gvs_06_mc_sim_ar1_blocks
```

---

## Interpretation

- FDR is expected to stay low (often well below the $\mathrm{tFDR}=0.1$ target) across the SNR grid for all three methods, but **TPR should recover much more slowly with SNR than in the equicorrelated Hastie/Mixed-Blocks demos** — since AR(1) decay means only nearby columns within a block are strongly correlated, weakening the group-level evidence aggregation that HAC clustering relies on.
- **IEN** is expected to sit between the near-zero TPR seen in the equicorrelated Gaussian demos and the EN-based methods: it should show real (though smaller) power gains as SNR increases, rather than staying flat.
- The $\rho$ sweep is the most informative view of this demo: watch how quickly recovery degrades as $\rho$ decreases and the AR(1) decay becomes steeper (correlation vanishing faster with column distance).

---

**Last updated**: 2026-07-08
