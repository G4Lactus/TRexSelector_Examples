# Demo 01: Screen-TRex on an i.i.d. Gaussian Design

## Purpose

Establish the **baseline behaviour of Screen-TRex** on the easiest possible design: independent Gaussian
 predictors with no correlation structure at all.
 Two thresholding rules are compared — the **Ordinary** majority-vote rule ($\Phi_j > 0.5$) and the
 **Bootstrap-CI** rule — over a signal-to-noise sweep.
 Everything the later demos add (correlation, dependency-aware corrections, biobank routing, solver
 backends) is measured against the numbers established here.
 Screening returns a *candidate set*, and FDR/TPR are evaluated on the individual selected variables
 (see [What is actually measured](../README.md#what-is-actually-measured-in-these-demos)).

---

## Data Generation Parameters (`make_iid_dgp`)

We consider the linear model:

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon},
\qquad \boldsymbol{\epsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma_{\varepsilon}^2 \boldsymbol{I}_n)
$$

- $\boldsymbol{y} \in \mathbb{R}^n$ is the response vector.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ is the design matrix.
- $\boldsymbol{\beta} \in \mathbb{R}^p$ is the coefficient vector, with $s$ nonzero entries.
- $\boldsymbol{\epsilon}$ is the noise vector, i.i.d. standard normal.
- $\sigma_{\varepsilon}^2$ is the noise variance, calibrated to achieve a target linear signal-to-noise ratio (SNR).
- $n = 300$, $p = 1000$, $s = 10$ (high-dimensional, $p > n$).

The design matrix has **no correlation structure**:

$$
X_{ij} \sim \mathcal{N}(0,1) \quad \text{i.i.d.}
$$

- The active support is drawn uniformly at random *per Monte Carlo trial*, so the results are not tied to
   one support pattern.
- All active coefficients are $\beta_j = 1$; the support-selection and coefficient RNGs are offset from the
   trial seed so they stay independent of the design and noise draws.

---

## Control Parameters

```text
K = 20                       # Random experiments per T-loop iteration
R_boot = 1000                # Bootstrap replicates (Bootstrap-CI rule only)
ci_grid_step = 0.001         # Bootstrap-CI threshold grid granularity
solver = TLARS               # T-Rex solver backend
MC = 200                     # Monte Carlo repetitions per grid point
```

Note that Screen-TRex has **no target-FDR parameter**: unlike the classical T-Rex selector, screening
 thresholds the voting statistic instead of calibrating to a user-specified level.

---

## Methods Compared

Two Screen-TRex thresholding rules [[1]](#references), both using `ScreenTRexMethod::TREX`:

- **Screen-TRex Ordinary** — selects $\{ j : \Phi_j > 0.5 \}$, a simple majority vote of the random
   experiments.
- **Screen-TRex Bootstrap** — builds a bootstrap confidence band around the estimated FDR curve
   (`R_boot = 1000` replicates) and picks its threshold from that band.

Both report an **estimated FDR** alongside the realized FDR/TPR — the procedure's own internal assessment
 of how many of its selections are false. Comparing the two is a central purpose of this suite.

---

## The Sweep

A single **SNR sweep** over $\mathrm{SNR} \in \{0.01, 0.1, 0.2, 0.5, 0.6, 1, 2, 5\}$, 200 MC trials per
point, with a fresh design, support, and noise draw in every trial.

---

## Output Files

Written to `simulation_results/data/`:

- `scr_screen_trex_snr_n300_p1000.txt` / `.csv` — FDR, TPR, and estimated FDR per method and SNR level.

Figures (PNG + PDF) go to `simulation_results/plots/`, produced by `./generate_plots.sh`.

---

## Running the Demo

```bash
./build/release/bin/trex_selector_methods/trex_screening/demo_trex_scr_01_mc_sim_screen_trex/demo_trex_scr_01_mc_sim_screen_trex
./generate_plots.sh   # render the figure below from the saved CSV
```

---

## Simulation Results

- **Screening is not FDR-controlling at low SNR.** At $\mathrm{SNR} \le 0.2$ the realized FDR runs between
   $0.16$ and $0.43$ for both rules. With almost no recoverable signal (TPR $\le 0.12$), nearly every
   selection is a false one — a property of thresholding a voting statistic without a target level, not a
   defect of the implementation.
- **Both rules become reliable once signal is present.** From $\mathrm{SNR} \ge 0.5$ the Ordinary rule's
   FDR drops to $0.12$ and settles around $0.05$–$0.07$, while its TPR climbs $0.44 \to 1.00$. All ten
   active variables are recovered at $\mathrm{SNR} = 5$.
- **Bootstrap-CI is the conservative rule: lower FDR, clearly less power.** Its FDR stays below $0.06$ for
   $\mathrm{SNR} \ge 0.5$ (roughly half the Ordinary rule's), paid for with TPR $0.87$ vs. $1.00$ at
   $\mathrm{SNR} = 5$, and a much wider gap mid-sweep ($0.37$ vs. $0.83$ at $\mathrm{SNR} = 1$). Choose it
   when false discoveries are costlier than missed ones.
- **The internal FDR estimate errs on the safe side here.** For $\mathrm{SNR} \ge 0.5$ it sits *above* the
   realized FDR for both rules (e.g. $0.115$ vs. $0.067$ for Ordinary at $\mathrm{SNR} = 1$) — it
   over-states the error rate rather than hiding it. At the lowest SNR points it becomes unreliable in both
   directions. Demo 03 shows this same estimate failing far more dramatically under correlation.

TPR (left) and FDR (right) vs. SNR (log axis), one line per thresholding rule; on the FDR panel the solid
line is the realized FDR and the dashed line the procedure's own estimated FDR.

![Screen-TRex on an i.i.d. design: TPR/FDR vs SNR](simulation_results/plots/scr_screen_trex_snr_n300_p1000.png)

---

## References

1. Machkour, J., Muma, M., & Palomar, D. P., "False Discovery Rate Control for Fast Screening of
   Large-Scale Genomics Biobanks.", IEEE Statistical Signal Processing Workshop (SSP), 2023,
    pp. 666–670, IEEE.
    [DOI-Link](https://doi.org/10.1109/SSP53291.2023.10207957)

---

**Last updated**: 2026-07-20
