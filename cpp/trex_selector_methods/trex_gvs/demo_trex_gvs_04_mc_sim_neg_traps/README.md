# Demo 04: T-Rex+GVS on Negative-Correlation Traps

## Purpose

Test whether T-Rex+GVS can correctly recover an **active group with negative within-group correlation structure** while **excluding two spatially separated inactive trap groups** that share the same sign-flipped correlation pattern. This stresses the grouping logic's ability to distinguish signal from confounding correlation geometry, rather than just correlation strength.

---

## Data Generation Parameters (`make_neg_traps_dgp`)

Each correlated group is built from a shared latent factor, with one half of the group loading **positively** and the other half loading **negatively** on it:

- **Active group $G_1$** (columns 0–99): positive half (0–49) $X_j = Z_1 + \sigma_x \xi_j$ with $\beta_j = +3$; negative half (50–99) $X_j = -Z_1 + \sigma_x \xi_j$ with $\beta_j = -3$.
- **Trap 1** (columns 100–199): same sign-flipped structure as $G_1$, but $\beta_j = 0$ (inactive).
- **Noise 1** (columns 200–299): i.i.d. $\mathcal{N}(0,1)$ background.
- **Trap 2** (columns 300–359): a second, smaller (60-column) sign-flipped trap, $\beta_j = 0$.
- **Noise 2** (columns 360–499): i.i.d. background.
- $n = 200$, $p = 500$, $s = 100$.
- Within-group correlation follows the same $\rho = 1/(1+\sigma_x^2)$ convention as the other block DGPs.

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

1. **SNR sweep** at fixed $\sigma_x = \sqrt{0.01}$: $\mathrm{SNR} \in \{0.1, 0.2, 0.5, 1.0, 2.0, 5.0\}$.
2. **$\rho$ sweep** at fixed $\mathrm{SNR}=2.0$ ($\sigma_x$ derived from $\rho$).
3. **2-D SNR $\times$ $\rho$ grid**.

---

## Output Files

Written to `simulation_results/`:

- `gvs_Neg-Traps_snr.txt` / `.csv`
- `gvs_Neg-Traps_rho.txt` / `.csv`
- `gvs_Neg-Traps_2d.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_gvs_04_mc_sim_neg_traps
```

---

## Interpretation

- **EN** and **EN+AUG** are expected to track each other closely, with FDR near the $\mathrm{tFDR}=0.1$ target and TPR reaching near-complete recovery already at moderate SNR — the sign-flipped active group is still identified correctly as a single group via the (unsigned) correlation-based HAC clustering step, and both spatially separated traps should be excluded.
- **IEN** is expected to remain conservative on this design (low FDR, low TPR that increases only slowly with SNR), similar to the uniformly-positive-correlation scenarios in Demos 01–03 — the sign-flipped halves here do not by themselves trigger substantially different IEN behavior.
- The two spatially separated traps are a stronger test of false-group exclusion than a single trap: watch whether FDR creeps up with SNR, which would indicate one of the traps starting to be (partially) selected.

---

**Last updated**: 2026-07-04
