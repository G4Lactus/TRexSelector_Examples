# Equi-Correlated Data

## Data Setup - Equi-Correlation

In this scenario, we generate $n$ observations on $p$ predictors through a
**single latent factor model**:

$$
X_{i,j} = \sqrt{\rho}\, f_i + \sqrt{1-\rho}\, \eta_{i,j},
$$
with $f_i \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0,1)$, and
$\eta_{i,j} \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0,1)$, and with all variables being mutually
independent.
Consequently, each row $\boldsymbol{x}_i = (X_{i,1},\ldots,X_{i,p})^\top$ is an $\text{i.i.d.}$ draw
from a zero-mean Gaussian with **compound symmetry** covariance

$$
\boldsymbol{\Sigma}_{\text{equi}} =
(1 - \rho) \, \boldsymbol{I}_p + \rho \, \mathbf{1}_p \mathbf{1}_p^{\top},
\quad \rho \in \left(-\frac{1}{p-1},\, 1\right).
$$

The diagonal entries satisfy 
$[\boldsymbol{\Sigma}_{\text{equi}}]_{jj}=1$ for every $j$, and the lower bound on
$\rho$ guarantees positive definiteness.
The eigenstructure is fully explicit:

$$
\begin{align*}
\lambda_1 &= 1 + (p-1)\rho && \text{with multiplicity 1 (eigenvector } \mathbf{1}_p/\sqrt{p}) \\
\lambda_2 &= 1 - \rho && \text{with multiplicity } (p-1)
\end{align*}
$$

When $\rho \leq \rho_{\text{thr\_DA}}$, the single factor is negligible.
As $\rho \to 1$, the $p - 1$ minor eigenvalues collapse to zero, the predictor matrix becomes
rank-1, and the regression is unidentifiable.

### Sparse Linear Model and SNR Control

Responses are generated from the standard **sparse linear model**

$$
\boldsymbol{y} =
\boldsymbol{X} \boldsymbol{\beta} + \boldsymbol{\varepsilon},
\qquad \boldsymbol{\varepsilon} \sim \mathcal{N}(0,\,\sigma^2 \boldsymbol{I}_n),
\quad \boldsymbol{\varepsilon} \perp \boldsymbol{X},
$$
$n$ as the number of observations, $p$ as the number of variables, and with a
**sparse coefficient vector**
$$
\boldsymbol{\beta} \in \mathbb{R}^p.
$$
Its support and sparsity level are

$$
\mathcal{S} = \operatorname{supp}(\beta) = \{j : \beta_j \neq 0\}, \qquad s = |\mathcal{S}| \ll p.
$$

The signal-to-noise ratio is defined based on
$$
\mathrm{SNR} =
\frac{\text{Var}(\boldsymbol{X} \boldsymbol{\beta})}{\sigma^{2}_{\varepsilon}}.
$$

## Simulation Setups

## Equi-Correlated SNR Sweep

We examine the scenario with $n=150$, $p=500$, $s=5$ actives,
$\boldsymbol{\beta}_{\text{active}} = 3$,
the support is generated according to the `Random` policy with a full redraw per Monte Carlo trial,
$\rho = \rho_{\text{between}} = \rho_{\text{within}} = 0.25$,
and $200$ Monte Carlo trials per SNR constellation, whereas the sweep of
$\text{SNR} \in \{0.10, 0.20, 0.50, 1.00, 2.00, 5.00\}$.

| SNR  | mean_FDP | mean_TPP | sd_FDP  | sd_TPP |
|----- |--------- |--------- |-------- |------- |
| 0.10 | 0.0000   | 0.0000   | 0.0000  | 0.0000 |
| 0.20 | 0.0000   | 0.0000   | 0.0000  | 0.0000 |
| 0.50 | 0.0000   | 0.0000   | 0.0000  | 0.0000 |
| 1.00 | 0.0000   | 0.0000   | 0.0000  | 0.0000 |
| 2.00 | 0.0000   | 0.0000   | 0.0000  | 0.0000 |
| 5.00 | 0.0000   | 0.0000   | 0.0000  | 0.0000 |

---

## Equi-Correlated $\rho$ Sweep

Keeping SNR constant while varying $n$, $p$, or $\rho$ isolates the statistical effect of the
equi-correlation structure from signal strength.

We examine the scenario with $n=150$, $p=500$, $s=5$ actives,
$\boldsymbol{\beta}_{\text{active}} = 3$,
the support is generated according to the `Random` policy with a full redraw per Monte Carlo trial.
The SNR is kept constant at $\text{SNR} = 2$, and the sweep takes place in
$\rho \in \{0.0, 0.1, 0.2, \ldots, 0.9 \}$.
As before $\rho = \rho_{\text{within}} = \rho_{\text{between}}$.

Results — 1D Pure Equi Sweep:

| rho  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|------|----------|----------|--------|------- |
| 0.00 | 0.1281   | 0.9955   | 0.1060 | 0.0208 |
| 0.10 | 0.0000   | 0.0000   | 0.0000 | 0.0000 |
| 0.20 | 0.0000   | 0.0000   | 0.0000 | 0.0000 |
| 0.30 | 0.0000   | 0.0000   | 0.0000 | 0.0000 |
| 0.40 | 0.0000   | 0.0000   | 0.0000 | 0.0000 |
| 0.50 | 0.0000   | 0.0000   | 0.0000 | 0.0000 |
| 0.60 | 0.0000   | 0.0000   | 0.0000 | 0.0000 |
| 0.70 | 0.0000   | 0.0000   | 0.0000 | 0.0000 |
| 0.80 | 0.0000   | 0.0000   | 0.0000 | 0.0000 |
| 0.90 | 0.0000   | 0.0000   | 0.0000 | 0.0000 |
