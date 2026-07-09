# `TREX+DA+AR(1)`

In the following we investigate the `TREX+DA+AR(1)` selector according to the version provided in
the CRAN package `TRexSelector`.

## Data Setup: AR(1)

The data generating process for an all-Toeplitz covariance based design matrix $\boldsymbol{X}$
assumes a **stationary first-order autoregressive process** along the predictor index.

$$
\Sigma_\text{AR1} = \begin{pmatrix} 1 & \rho & \rho^2 & \cdots & \rho^{p-1} \\ \rho & 1 & \rho & \cdots & \rho^{p-2} \\ \vdots & & \ddots & & \vdots \\ \rho^{p-1} & \cdots & & \rho & 1 \end{pmatrix}, \quad [\Sigma_\text{AR1}]_{jk} = \rho^{|j-k|} \, .
$$

The rows of $\boldsymbol{X}$ are generated as:

$$
X_{i,j} = \rho\, X_{i,j-1} + \eta_{i,j}, \quad \eta_{i,j} \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0,\, 1-\rho^2), \quad |\rho| < 1
$$
initialised at $X_{i,0} \sim \mathcal{N}(0,1)$.
The spectral density of this process is:
$$
S(\omega) = \frac{1 - \rho^2}{1 - 2\rho\cos\omega + \rho^2}, \quad \omega \in [-\pi, \pi]
$$
which has a single peak at $\omega = 0$ for $\rho > 0$ (low-frequency dominance) and at $\omega = \pm\pi$ for $\rho < 0$ (high-frequency dominance). The precision matrix is tridiagonal:

$$
\Omega_\text{AR1} = \frac{1}{1-\rho^2} \begin{pmatrix} 1 & -\rho & 0 & \cdots \\ -\rho & 1+\rho^2 & -\rho & \\ 0 & -\rho & \ddots & \\ \vdots & & & 1 \end{pmatrix} \, .
$$

The window half-width $\kappa = \lceil \log \rho_\text{thr} / \log |\rho| \rceil$ is the smallest
lag at which $|\rho|^\kappa \leq \rho_\text{thr}$, truncating the infinite-range correlation to
a finite sliding window.

## Simulation Setups

## SNR Sweep Scenario

We examine the performance of the selector for the SNR sweep
$\text{SNR} \in \{0.10, 0.20, ..., 1, 2, 5\}$.
The target FDR $\alpha = \text{tFDR} = 0.20$, $n=300$, $p=1000$, $s=10$,
$\boldsymbol{\beta}_{\text{active}}=3$, $\rho=0.7$, and $200$ Monte Carlo trials (MCs).

### Random Support per trial

The first scenario considers a random support for the active variables which changes per Monte
Carlo trials.

| SNR  | mean FDP | mean TPP | sd FDP | sd TPP |
|------|----------|----------|--------|--------|
| 0.10 | 0.0175   | 0.0020   | 0.1266 | 0.0172 |
| 0.20 | 0.0611   | 0.0225   | 0.2072 | 0.0525 |
| 0.50 | 0.0971   | 0.1955   | 0.2034 | 0.1717 |
| 1.00 | 0.0574   | 0.4360   | 0.1270 | 0.2157 |
| 2.00 | 0.0781   | 0.6055   | 0.1301 | 0.1931 |
| 5.00 | 0.1054   | 0.6205   | 0.1373 | 0.1906 |

### CappedSpread support

The second scenario assumes CappedSpread support with a maximum gap `max_gap=20`.
The support indices are then generated as $\{1, 21, 41, 61, ...\}$.

| SNR  | mean FDP | mean TPP | sd FDP | sd TPP |
|------|----------|----------|--------|--------|
| 0.10 | 0.0242   | 0.0060   | 0.1460 | 0.0277 |
| 0.20 | 0.0599   | 0.0250   | 0.2026 | 0.0616 |
| 0.50 | 0.1198   | 0.2215   | 0.1986 | 0.1629 |
| 1.00 | 0.0680   | 0.5660   | 0.1106 | 0.1703 |
| 2.00 | 0.0587   | 0.7910   | 0.0891 | 0.1563 |
| 5.00 | 0.1071   | 0.7935   | 0.1201 | 0.1632 |

## Rho Sweep

We examine the performance of the selector for the $\rho$ sweep
$\rho \in \{0.1, 0.2, \ldots, 0.9\}$.
The target FDR $\alpha := \text{tFDR} = 0.20$, $n=300$, $p=1000$, $s=10$, and $\text{SNR} = 2$.
The support is chosen according to the `CappedSpread` policy with an offset of $10$ over all
$200$ Monte Carlo trials.

| $\rho$ | mean FDP | mean TPP | sd FDP | sd TPP |
| ------ | ---------- | ---------- | -------- | -------- |
| 0.10 | 0.0866 | 0.9830 | 0.0910 | 0.0415 |
| 0.20 | 0.0778 | 0.9755 | 0.0858 | 0.0526 |
| 0.30 | 0.0734 | 0.9650 | 0.0812 | 0.0608 |
| 0.40 | 0.0788 | 0.9410 | 0.0863 | 0.0834 |
| 0.50 | 0.0816 | 0.9075 | 0.0951 | 0.1065 |
| 0.60 | 0.0676 | 0.8650 | 0.0919 | 0.1223 |
| 0.70 | 0.0587 | 0.7910 | 0.0891 | 0.1563 |
| 0.80 | 0.0624 | 0.5395 | 0.1303 | 0.1912 |
| 0.90 | 0.0950 | 0.0000 | 0.2940 | 0.0000 |

## SNR vs. Rho sweep

### The boundary case: `max_gap=1`

With `max_gap = 1`, the support is $\{1, 2, 3, …, 10\}$ and so consecutive active variables are
separated by only one index step.
Thus their pairwise correlation is $\rho^{1} = \rho$ itself.
That's the strongest possible inter-active correlation under AR(1).
This scenario results into an erosion such that the DA correction computes neighbourhood windows
from the estimated $\rho$ and deflates relative occurrences within those windows.
When active variables are inside each other's windows, the deflation suppresses their own votes,
making it harder for them to accumulate the relative occurrence mass required to cross the voting
threshold.
Consequently, the selector then either selects:

* None of them (TPP → 0 at high rho) or
* selects the wrong neighbour (FP)

This makes the max_gap parameter effectively a signal separation axis that is orthogonal to SNR.

Thus we examine a combination of a gap sweep with a rho sweep.
For each rho there is a minimum gap below which the DA correction can no longer recover the full
support.

The `TREX+DA+AR1` correction computes a window half-width:
$$
\kappa = ⌈ \log(\rho_{\text{thr}}) /  \log(\rho)⌉, \quad \rho_{\text{thr}} = 0.02 \text{(default)}.
$$
Active variables with $\text{gap} < \kappa$ fall inside each other's correction windows —
deflation fires on true positives and TPP collapses. The critical boundary is $\text{gap} = \kappa$:

| ρ   | κ  | gap=100 | gap=50 | gap=20   | gap=15   | gap=10   | gap=5    | gap=1 |
|-----|----|---------|--------|----------|----------|----------|----------|-------|
| 0.1 |  2 | ✓       | ✓      | ✓        | ✓        | ✓        | ✓        | ✗     |
| 0.2 |  3 | ✓       | ✓      | ✓        | ✓        | ✓        | ✓        | ✗     |
| 0.3 |  4 | ✓       | ✓      | ✓        | ✓        | ✓        | ✓        | ✗     |
| 0.4 |  5 | ✓       | ✓      | ✓        | ✓        | ✓        | boundary | ✗     |
| 0.5 |  6 | ✓       | ✓      | ✓        | ✓        | ✓        | ✗        | ✗     |
| 0.6 |  8 | ✓       | ✓      | ✓        | ✓        | ✗        | ✗        | ✗     |
| 0.7 | 11 | ✓       | ✓      | ✓        | ✗        | ✗        | ✗        | ✗     |
| 0.8 | 18 | ✓       | ✓      | ✗        | ✗        | ✗        | ✗        | ✗     |
| 0.9 | 38 | ✓       | ✓      | ✗        | ✗        | ✗        | ✗        | ✗     |

### 2D Gap vs. Rho Sweep

In the first scenario we consider $\alpha := \text{tFDR} = 0.20$, $n=300$, $p=1000$, $s=10$,
$SNR=2.0$, $200 MC$, and we evluate the gap grid `gap_grid  : {100, 50, 20, 15, 10, 5, 1}` and the
rho grid `rho_grid  : {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}`.

### DA+AR(1) Correction Window

κ = ⌈log(0.02) / log(ρ)⌉

| ρ     | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 |
|-------|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| **κ** |  2  |  3  |  4  |  5  |  6  |  8  | 11  | 18  | 38  |

---

## Scenario $\text{tFDR}=0.2$: mean TPP (CappedSpread gap values + Random)

| ρ       | gap=100 | gap=50 | gap=20 | gap=15 | gap=10 | gap=5  | gap=1  | Random |
|---------|---------|--------|--------|--------|--------|--------|--------|--------|
| **0.1** | 0.981   | 0.982  | 0.983  | 0.982  | 0.983  | 0.986  | 0.000  | 0.959  |
| **0.2** | 0.973   | 0.975  | 0.976  | 0.972  | 0.975  | 0.978  | 0.000  | 0.921  |
| **0.3** | 0.958   | 0.962  | 0.965  | 0.962  | 0.959  | 0.954  | 0.000  | 0.888  |
| **0.4** | 0.939   | 0.942  | 0.941  | 0.935  | 0.938  | 0.000  | 0.000  | 0.849  |
| **0.5** | 0.914   | 0.905  | 0.907  | 0.896  | 0.907  | 0.000  | 0.000  | 0.794  |
| **0.6** | 0.853   | 0.854  | 0.865  | 0.849  | 0.824  | 0.000  | 0.000  | 0.715  |
| **0.7** | 0.759   | 0.782  | 0.791  | 0.739  | 0.000  | 0.000  | 0.000  | 0.606  |
| **0.8** | 0.645   | 0.650  | 0.539  | 0.000  | 0.000  | 0.000  | 0.000  | 0.454  |
| **0.9** | 0.448   | 0.436  | 0.000  | 0.000  | 0.000  | 0.000  | 0.000  | 0.217  |

---

## Scenario $\text{tFDR}=0.2$: mean FDP (CappedSpread gap values + Random)

| ρ       | gap=100 | gap=50 | gap=20 | gap=15 | gap=10 | gap=5  | gap=1  | Random |
|---------|---------|--------|--------|--------|--------|--------|--------|--------|
| **0.1** | 0.086   | 0.086  | 0.087  | 0.076  | 0.083  | 0.089  | 0.025  | 0.083  |
| **0.2** | 0.085   | 0.087  | 0.078  | 0.078  | 0.080  | 0.090  | 0.025  | 0.078  |
| **0.3** | 0.083   | 0.082  | 0.073  | 0.080  | 0.081  | 0.092  | 0.020  | 0.080  |
| **0.4** | 0.080   | 0.084  | 0.079  | 0.080  | 0.080  | 0.005  | 0.020  | 0.073  |
| **0.5** | 0.081   | 0.089  | 0.082  | 0.079  | 0.085  | 0.015  | 0.025  | 0.077  |
| **0.6** | 0.078   | 0.090  | 0.068  | 0.066  | 0.082  | 0.010  | 0.020  | 0.087  |
| **0.7** | 0.079   | 0.080  | 0.059  | 0.061  | 0.015  | 0.020  | 0.020  | 0.078  |
| **0.8** | 0.062   | 0.079  | 0.062  | 0.095  | 0.050  | 0.055  | 0.065  | 0.066  |
| **0.9** | 0.095   | 0.107  | 0.095  | 0.095  | 0.110  | 0.085  | 0.100  | 0.088  |

---

The second scenario is almost identical, but $\alpha := \text{tFDR} = 0.10$.

## Scenario $\text{tFDR}=0.1$: mean TPP (CappedSpread gap values + Random)

| ρ       | gap=100 | gap=50 | gap=20 | gap=15 | gap=10 | gap=5  | gap=1  | Random |
|---------|---------|--------|--------|--------|--------|--------|--------|--------|
| **0.0** | 0.950   | 0.965  | 0.967  | 0.957  | 0.963  | 0.967  | 0.001  | 0.960  |
| **0.1** | 0.953   | 0.961  | 0.967  | 0.956  | 0.964  | 0.965  | 0.000  | 0.936  |
| **0.2** | 0.946   | 0.958  | 0.959  | 0.943  | 0.959  | 0.960  | 0.000  | 0.900  |
| **0.3** | 0.934   | 0.939  | 0.943  | 0.937  | 0.946  | 0.940  | 0.000  | 0.866  |
| **0.4** | 0.918   | 0.930  | 0.922  | 0.906  | 0.927  | 0.000  | 0.000  | 0.826  |
| **0.5** | 0.890   | 0.889  | 0.887  | 0.885  | 0.881  | 0.000  | 0.000  | 0.788  |
| **0.6** | 0.835   | 0.838  | 0.833  | 0.824  | 0.822  | 0.000  | 0.000  | 0.711  |
| **0.7** | 0.742   | 0.763  | 0.760  | 0.730  | 0.000  | 0.000  | 0.000  | 0.594  |
| **0.8** | 0.625   | 0.669  | 0.551  | 0.000  | 0.000  | 0.000  | 0.000  | 0.438  |
| **0.9** | 0.434   | 0.424  | 0.000  | 0.000  | 0.000  | 0.000  | 0.000  | 0.201  |

---

## Scenario $\text{tFDR}=0.1$: mean FDP (CappedSpread gap values + Random)

| ρ       | gap=100 | gap=50 | gap=20 | gap=15 | gap=10 | gap=5  | gap=1  | Random |
|---------|---------|--------|--------|--------|--------|--------|--------|--------|
| **0.0** | 0.027   | 0.027  | 0.024  | 0.026  | 0.017  | 0.029  | 0.015  | 0.031  |
| **0.1** | 0.026   | 0.028  | 0.027  | 0.024  | 0.024  | 0.026  | 0.010  | 0.021  |
| **0.2** | 0.030   | 0.027  | 0.028  | 0.028  | 0.027  | 0.030  | 0.015  | 0.026  |
| **0.3** | 0.030   | 0.029  | 0.031  | 0.033  | 0.026  | 0.026  | 0.015  | 0.030  |
| **0.4** | 0.029   | 0.030  | 0.026  | 0.033  | 0.024  | 0.005  | 0.010  | 0.034  |
| **0.5** | 0.037   | 0.035  | 0.025  | 0.026  | 0.034  | 0.005  | 0.015  | 0.028  |
| **0.6** | 0.027   | 0.042  | 0.032  | 0.028  | 0.027  | 0.005  | 0.010  | 0.030  |
| **0.7** | 0.031   | 0.038  | 0.024  | 0.033  | 0.015  | 0.020  | 0.015  | 0.025  |
| **0.8** | 0.033   | 0.023  | 0.034  | 0.090  | 0.050  | 0.055  | 0.065  | 0.026  |
| **0.9** | 0.096   | 0.094  | 0.085  | 0.095  | 0.105  | 0.080  | 0.095  | 0.078  |
