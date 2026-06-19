# Simulation Study: `TREX+DA+BT` with AR(1) Block Correlated Data

## Setup

The following simulations examine the performance of the `TREX+DA+BT` selector for block-diagonal
data consisting of $M$ AR(1) blocks depending on a homogeneous autocorrelation $\rho$ for all blocks
within the linear model
$$
\boldsymbol{y} = \boldsymbol{X} \boldsymbol{\beta} + \boldsymbol{\varepsilon},
$$
with $n$ observations and $p = M \cdot Q = p_{\text{ar}}$ variables, i.e., the data are split into 
$M$ groups of equal size $Q$.
The simulations investigate ceteris-paribus (as far as applicable) effects of varying one
parameter while fixing all others across **low- to high-dimensional regimes**.

The **within**-group correlation structure of the AR(1) blocks follows the structure of the
$Q \times Q$ Toeplitz matrix
$$
\boldsymbol{\Sigma_{m}}(\rho) =
\begin{bmatrix}
1 & \rho & \rho^{2} & \cdots & \rho^{Q-1} \\
\rho & 1 & \rho & \cdots & \rho^{Q-2} \\
\rho^{2} & \rho & 1 & \ldots & \rho^{Q-3} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
\rho^{Q-1} & \rho^{Q-2} & \rho^{Q-3} & \cdots & 1
\end{bmatrix} \in \mathbb{R}^{Q \times Q}.
$$
with
$$
[\boldsymbol{\Sigma}_{m}(\rho)]_{j,k} = \rho^{|j - k|}.
$$

The data generation **proceeds via** recursion for each block $\boldsymbol{X}^{(m)}_i$
$$
z_{j} = \rho \, z_{j-1} + \sqrt{1 - \rho^{2}} \, \varepsilon_{j}, \quad j=2,\ldots,Q.
$$

The **total covariance** is then
$$
\boldsymbol{\Sigma} =
\begin{bmatrix}
\boldsymbol{\Sigma}_1(\rho) & \mathbf{0} & \cdots & \mathbf{0}\\
\mathbf{0} & \boldsymbol{\Sigma}_2(\rho) & \cdots & \mathbf{0}\\
\vdots & \vdots & \ddots & \vdots\\
\mathbf{0} & \mathbf{0} & \cdots & \boldsymbol{\Sigma}_M(\rho)
\end{bmatrix} \in \mathbb{R}^{p \times p}.
$$

The $i$-th row of $\boldsymbol{X}$ is denoted as column vector with
$\boldsymbol{X}_i^{(m)} \in \mathbb{R}^Q$ comprising blocks
$$
\boldsymbol{X}_i = \begin{pmatrix}
\boldsymbol{X}_i^{(1)} \\ \boldsymbol{X}_i^{(2)} \\ \vdots \\
\boldsymbol{X}_i^{(M)}
\end{pmatrix}.
$$

### Support Considerations

The true support $S$ draws **exactly one representative per AR(1) block**:
$$
S = \{s_1, \dots, s_M\}, \quad s_m \in \{(m-1)Q + 1, \dots, mQ\},
$$
such that the active indices within $\boldsymbol{\beta}$ are
$$
\beta_{s_m} = \beta_{\text{active}} > 0, \, \forall m = 1,\ldots,M.
$$

---

## 1. AR(1) Blocks (no white block): SNR Sweep

We consider $n = 150$, $M = 5$, $Q = 5$, $p = M * Q = 25$, $s=5$, $\rho = 0.7$, $200$ MCs,
$\text{tFDR} = 0.20$ and $\text{SNR} \in \{0.1, 0.2, 0.5, 0.6, 1, 2, 5\}$,
$\beta_{\text{active}} = 1$.

### 1a.: `hc_dist = "single"`

|  SNR  | mean_FDP | mean_TPP | sd_FDP  | sd_TPP |
|-------|----------|----------|---------|--------|
|  0.10 | 0.0363   | 0.0290   | 0.1779  | 0.1068 |
|  0.20 | 0.0604   | 0.0990   | 0.2043  | 0.2015 |
|  0.50 | 0.0851   | 0.4330   | 0.1811  | 0.3437 |
|  0.60 | 0.1026   | 0.5350   | 0.1996  | 0.3443 |
|  1.00 | 0.0877   | 0.7270   | 0.1575  | 0.2968 |
|  2.00 | 0.0847   | 0.7750   | 0.1532  | 0.2926 |
|  5.00 | 0.0852   | 0.7760   | 0.1523  | 0.2916 |

### 1b.: `hc_dist = "complete"`

| SNR  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|------|----------|----------|--------|--------|
| 0.10 | 0.0488   | 0.0360   | 0.2034 | 0.1199 |
| 0.20 | 0.0901   | 0.0990   | 0.2405 | 0.1975 |
| 0.50 | 0.1117   | 0.5050   | 0.1970 | 0.3426 |
| 0.60 | 0.1300   | 0.5860   | 0.1875 | 0.3328 |
| 1.00 | 0.1377   | 0.7740   | 0.1744 | 0.2507 |
| 2.00 | 0.1416   | 0.8130   | 0.1675 | 0.2296 |
| 5.00 | 0.1460   | 0.8110   | 0.1721 | 0.2331 |

### 1c.: `hc_dist = "average"`

| SNR  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|------|----------|----------|--------|--------|
| 0.10 | 0.0538   | 0.0350   | 0.2142 | 0.1194 |
| 0.20 | 0.0926   | 0.0980   | 0.2427 | 0.1933 |
| 0.50 | 0.1079   | 0.5030   | 0.1895 | 0.3373 |
| 0.60 | 0.1276   | 0.5820   | 0.1873 | 0.3290 |
| 1.00 | 0.1358   | 0.7650   | 0.1783 | 0.2610 |
| 2.00 | 0.1417   | 0.8050   | 0.1780 | 0.2443 |
| 5.00 | 0.1440   | 0.8040   | 0.1771 | 0.2455 |

---

## 2. AR(1) Block (no white block): $\rho$ Sweep

We consider $n = 150$, $M = 5$, $Q = 5$, $p = M * Q = 25$, $s = 5$, $\text{SNR} = 2.0$, $200$ MCs,
$\text{tFDR} = 0.20$ and $\rho \in \{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9\}$.

### 2a.: `hc_dist = "single"`

| rho | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|-----|----------|----------|--------|--------|
| 0.0 | 0.0318   | 0.4090   | 0.1097 | 0.2578 |
| 0.1 | 0.0339   | 0.4000   | 0.1207 | 0.2630 |
| 0.2 | 0.0407   | 0.4790   | 0.1221 | 0.2926 |
| 0.3 | 0.0615   | 0.6830   | 0.1312 | 0.3043 |
| 0.4 | 0.0559   | 0.8520   | 0.1368 | 0.2383 |
| 0.5 | 0.0632   | 0.8580   | 0.1422 | 0.2229 |
| 0.6 | 0.0826   | 0.8170   | 0.1588 | 0.2762 |
| 0.7 | 0.0847   | 0.7750   | 0.1532 | 0.2926 |
| 0.8 | 0.1479   | 0.6800   | 0.2142 | 0.3080 |
| 0.9 | <span style="background-color: red; color: white;">0.2192</span>   | 0.5760   | 0.2360 | 0.3005 |

### 2a.: `hc_dist = "complete"`

| rho | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|-----|----------|----------|--------|--------|
| 0.0 | 0.0756   | 0.7350   | 0.1307 | 0.2422 |
| 0.1 | 0.0846   | 0.7250   | 0.1441 | 0.2817 |
| 0.2 | 0.0707   | 0.7620   | 0.1160 | 0.2803 |
| 0.3 | 0.0973   | 0.7760   | 0.1484 | 0.2664 |
| 0.4 | 0.0821   | 0.8160   | 0.1338 | 0.2483 |
| 0.5 | 0.0889   | 0.8640   | 0.1518 | 0.2120 |
| 0.6 | 0.1117   | 0.8480   | 0.1628 | 0.2235 |
| 0.7 | 0.1416   | 0.8130   | 0.1675 | 0.2296 |
| 0.8 | <span style="background-color: red; color: white;">0.2023</span>   | 0.7660   | 0.1883 | 0.2373 |
| 0.9 | <span style="background-color: red; color: white;">0.2652</span>   | 0.6870   | 0.2165 | 0.2549 |

### 2a.: `hc_dist = "average"`

| rho | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|-----|----------|----------|--------|------- |
| 0.0 | 0.0701   | 0.6810   | 0.1282 | 0.2529 |
| 0.1 | 0.0739   | 0.6990   | 0.1319 | 0.2836 |
| 0.2 | 0.0795   | 0.7190   | 0.1290 | 0.2864 |
| 0.3 | 0.0903   | 0.8170   | 0.1501 | 0.2387 |
| 0.4 | 0.0702   | 0.8790   | 0.1218 | 0.2014 |
| 0.5 | 0.0765   | 0.8850   | 0.1436 | 0.1989 |
| 0.6 | 0.1049   | 0.8480   | 0.1528 | 0.2253 |
| 0.7 | 0.1418   | 0.8050   | 0.1814 | 0.2459 |
| 0.8 | <span style="background-color: red; color: white;">0.2052</span>   | 0.7450   | 0.1930 | 0.2520 |
| 0.9 | <span style="background-color: red; color: white;">0.2570</span>   | 0.6790   | 0.2165 | 0.2567 |

---

## 3. AR(1) Block: Q sweep

We consider $n = 150$, $M = 5$, $s = 5$, $\rho = 0.7$, $\text{SNR} = 2.0$, $200$ MCs,
$Q \in \{5, 10, 15, 20, 25, 30, 35, 40, 45, 50\}$.
Thus, $p = M * Q $ varies from 25 to 250.

### 3a.: `hc_dist = "single"`

| Q  | p   | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|-----|----------|----------|--------|--------|
| 5  | 25  | 0.0847   | 0.7750   | 0.1532 | 0.2926 |
| 10 | 50  | 0.1184   | 0.6870   | 0.1858 | 0.3406 |
| 15 | 75  | 0.1190   | 0.6810   | 0.1724 | 0.3570 |
| 20 | 100 | 0.1250   | 0.6760   | 0.1716 | 0.3533 |
| 25 | 125 | 0.1405   | 0.6540   | 0.1986 | 0.3579 |
| 30 | 150 | 0.1322   | 0.6630   | 0.1647 | 0.3435 |
| 35 | 175 | 0.1380   | 0.6510   | 0.1916 | 0.3496 |
| 40 | 200 | 0.1400   | 0.6620   | 0.1794 | 0.3327 |
| 45 | 225 | 0.1337   | 0.6230   | 0.1810 | 0.3703 |
| 50 | 250 | 0.1230   | 0.6630   | 0.1816 | 0.3376 |

### 3b.: `hc_dist = "complete"`
  
| Q  | p   | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|-----|----------|----------|--------|--------|
| 5  | 25  | 0.1416   | 0.8130   | 0.1675 | 0.2296 |
| 10 | 50  | 0.1437   | 0.7800   | 0.1741 | 0.2919 |
| 15 | 75  | 0.1569   | 0.7900   | 0.1654 | 0.2762 |
| 20 | 100 | 0.1668   | 0.8100   | 0.1707 | 0.2510 |
| 25 | 125 | 0.1644   | 0.7760   | 0.1809 | 0.2923 |
| 30 | 150 | 0.1596   | 0.7740   | 0.1880 | 0.2942 |
| 35 | 175 | 0.1542   | 0.7940   | 0.1804 | 0.2704 |
| 40 | 200 | 0.1433   | 0.7910   | 0.1598 | 0.2751 |
| 45 | 225 | 0.1457   | 0.7970   | 0.1699 | 0.2723 |
| 50 | 250 | 0.1471   | 0.7690   | 0.1853 | 0.2822 |

### 3c.: `hc_dist = "average"`

| Q  | p   | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|-----|----------|----------|--------|--------|
| 5  | 25  | 0.1417   | 0.8050   | 0.1780 | 0.2443 |
| 10 | 50  | 0.1439   | 0.7780   | 0.1751 | 0.2911 |
| 15 | 75  | 0.1527   | 0.7920   | 0.1661 | 0.2792 |
| 20 | 100 | 0.1689   | 0.7990   | 0.1823 | 0.2641 |
| 25 | 125 | 0.1649   | 0.7720   | 0.1791 | 0.2906 |
| 30 | 150 | 0.1521   | 0.7790   | 0.1773 | 0.2881 |
| 35 | 175 | 0.1492   | 0.7940   | 0.1777 | 0.2689 |
| 40 | 200 | 0.1394   | 0.7900   | 0.1590 | 0.2827 |
| 45 | 225 | 0.1525   | 0.7960   | 0.1755 | 0.2690 |
| 50 | 250 | 0.1423   | 0.7710   | 0.1852 | 0.2852 |

---

## 4. AR(1) Block: M sweep

We consider $n=150$, $Q=5$, $\rho=0.7$, $\text{SNR}=2.0$, $200$ MCs,
$M \in \{1, 3, 5, 10, 15, 20, 30\}$, thus $p = M * Q$, s = M varies.

### 4a.: `hc_dist = "single"`

| M  | p   | s  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|-----|----|----------|----------|--------|--------|
| 1  | 5   | 1  | 0.0075   | 0.1850   | 0.0789 | 0.3893 |
| 3  | 15  | 3  | 0.1083   | 0.5000   | 0.2053 | 0.4227 |
| 5  | 25  | 5  | 0.0847   | 0.7750   | 0.1532 | 0.2926 |
| 10 | 50  | 10 | 0.0772   | 0.6645   | 0.1262 | 0.2235 |
| 15 | 75  | 15 | 0.1159   | 0.3457   | 0.1745 | 0.2153 |
| 20 | 100 | 20 | 0.0984   | 0.1477   | 0.1569 | 0.1490 |
| 30 | 150 | 30 | 0.0720   | 0.0405   | 0.1818 | 0.0725 |

### 4a.: `hc_dist = "complete"`

| M  | p   | s  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|-----|----|----------|----------|--------|--------|
| 1  | 5   | 1  | 0.0075   | 0.1850   | 0.0789 | 0.3893 |
| 3  | 15  | 3  | 0.1294   | 0.5783   | 0.1962 | 0.4031 |
| 5  | 25  | 5  | 0.1416   | 0.8130   | 0.1675 | 0.2296 |
| 10 | 50  | 10 | 0.1138   | 0.6805   | 0.1400 | 0.2395 |
| 15 | 75  | 15 | 0.1460   | 0.3970   | 0.1536 | 0.2302 |
| 20 | 100 | 20 | 0.1380   | 0.1880   | 0.1727 | 0.1628 |
| 30 | 150 | 30 | 0.1224   | 0.0547   | 0.2276 | 0.0833 |

### 4a.: `hc_dist = "average"`

| M  | p   | s  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|-----|----|----------|----------|--------|--------|
| 1  | 5   | 1  | 0.0075   | 0.1850   | 0.0789 | 0.3893 |
| 3  | 15  | 3  | 0.1199   | 0.5800   | 0.1866 | 0.4082 |
| 5  | 25  | 5  | 0.1417   | 0.8050   | 0.1780 | 0.2443 |
| 10 | 50  | 10 | 0.1129   | 0.6715   | 0.1427 | 0.2386 |
| 15 | 75  | 15 | 0.1459   | 0.3923   | 0.1551 | 0.2300 |
| 20 | 100 | 20 | 0.1351   | 0.1867   | 0.1692 | 0.1618 |
| 30 | 150 | 30 | 0.1165   | 0.0532   | 0.2261 | 0.0831 |
