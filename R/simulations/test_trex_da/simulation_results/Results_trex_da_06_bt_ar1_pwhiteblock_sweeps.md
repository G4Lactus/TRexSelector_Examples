# Simulation Study: `TREX+DA+BT` with AR(1) Block Correlated Data + White Block Data

## Setup

The following simulations examine the performance of the `TREX+DA+BT` selector for block-diagonal
data consisting of $M$ AR(1) blocks depending on a homogeneous autocorrelation $\rho$ for all blocks
and one block of white variables to compensate for a high-dimensional ($p>n$) variable selection
scenario within the linear model
$$
\boldsymbol{y} = \boldsymbol{X} \boldsymbol{\beta} + \boldsymbol{\varepsilon}.
$$
We consider $n$ observations and $p = p_{\text{AR}} + p_{\text{white}}$ variables.
The simulations investigate ceteris-paribus (as far as applicable) effects of varying one
parameter while fixing all others.

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
\end{bmatrix} \in \mathbb{R}^{Q \times Q},
$$
with
$$
[\boldsymbol{\Sigma}_{m}(\rho)]_{j,k} = \rho^{|j - k|}.
$$

In total there are $M$ statistically independent AR(1) groups with $Q$ members each, thus
$p_{\text{ar}} = M \cdot Q$.
The data generation is achieved for each AR(1) block $\boldsymbol{X}^{(m)}_{i}$ using the recursion:
$$
z_{j} = \rho \, z_{j-1} + \sqrt{1 - \rho^{2}} \, \varepsilon_{j},
$$
for $j = 2, \ldots, Q$ with $z_{1}, \varepsilon_{j} \sim \mathcal{N}(0, 1)$, for $m = 1, \ldots, M$.

The $M+1$ block represents white variables $\boldsymbol{X}^{(M+1)}_{i} \sim \mathcal{N}(0, 1)$
$\forall j' = 1, \ldots, p_{\text{white}}$, where $p_{\text{white}} = p - p_{\text{ar}}$.

The **total covariance** is then
$$
\boldsymbol{\Sigma} =
\begin{bmatrix}
\boldsymbol{\Sigma}_1(\rho) & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} \\
\mathbf{0} & \boldsymbol{\Sigma}_2(\rho) & \cdots & \mathbf{0} & \mathbf{0} \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
\mathbf{0} & \mathbf{0} & \cdots & \boldsymbol{\Sigma}_M(\rho) & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{I}_{p_\text{white}}
\end{bmatrix} \in \mathbb{R}^{p \times p}.
$$
The $i$-th row of $\boldsymbol{X}$ is denoted as column vector with
$\boldsymbol{X}_i^{(m)} \in \mathbb{R}^Q$ comprising
$$
\boldsymbol{X}_i = \begin{pmatrix}
\boldsymbol{X}_i^{(1)} \\ \boldsymbol{X}_i^{(2)} \\ \vdots \\ \boldsymbol{X}_i^{(M)} \\
\boldsymbol{X}_{i}^{(M+1)}
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

## 1. AR(1) Blocks + White Block: SNR sweep

We consider $n=150$, $p=1000$, $M=5$, $Q=5$, $p_{\text{ar}}=25$, $p_{\text{white}}=975$, $s=5$,
$\rho=0.7$, $200$ MCs, $\text{SNR} \in \{0.1, 0.2, 0.5, 0.6, 1, 2, 5\}$.

### 1a.: `hc_dist = "single"`

| SNR  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|------|----------|----------|--------|--------|
| 0.10 | 0.0000   | 0.0010   | 0.0000 | 0.0141 |
| 0.20 | 0.0025   | 0.0100   | 0.0354 | 0.0558 |
| 0.50 | 0.0300   | 0.1390   | 0.1423 | 0.1862 |
| 0.60 | 0.0331   | 0.2410   | 0.1069 | 0.2409 |
| 1.00 | 0.0324   | 0.5760   | 0.0925 | 0.2739 |
| 2.00 | 0.0561   | 0.8570   | 0.1370 | 0.1948 |
| 5.00 | 0.0733   | 0.8910   | 0.1599 | 0.1791 |

### 1b.: `hc_dist = "complete"`

| SNR  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|------|----------|----------|--------|--------|
| 0.10 | 0.0000   | 0.0010   | 0.0000 | 0.0141 |
| 0.20 | 0.0045   | 0.0100   | 0.0452 | 0.0558 |
| 0.50 | 0.0618   | 0.1630   | 0.1937 | 0.2006 |
| 0.60 | 0.0622   | 0.2680   | 0.1433 | 0.2484 |
| 1.00 | 0.0780   | 0.6090   | 0.1346 | 0.2692 |
| 2.00 | 0.1027   | 0.8650   | 0.1453 | 0.1628 |
| 5.00 | 0.1118   | 0.8940   | 0.1514 | 0.1552 |

### 1c.: `hc_dist = "average"`

| SNR  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|------|----------|----------|--------|--------|
| 0.10 | 0.0000   | 0.0010   | 0.0000 | 0.0141 |
| 0.20 | 0.0050   | 0.0090   | 0.0499 | 0.0461 |
| 0.50 | 0.0606   | 0.1570   | 0.1939 | 0.1953 |
| 0.60 | 0.0589   | 0.2620   | 0.1396 | 0.2434 |
| 1.00 | 0.0764   | 0.6030   | 0.1381 | 0.2701 |
| 2.00 | 0.1007   | 0.8550   | 0.1485 | 0.1724 |
| 5.00 | 0.1128   | 0.8840   | 0.1576 | 0.1658 |

---

## 2. AR(1) Blocks + White Block: $\rho$ sweep

We consider $n=150$, $p=1000$, $M=5$, $Q=5$, $p_ar=25$, $p_white=975$, $s=5$, $\text{SNR}=2.0$,
$200$ MCs, $\rho \in \{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9\}$.

### 2a.: `hc_dist = "single"`

| rho | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|-----|----------|----------|--------|--------|
| 0.0 | 0.0603   | 0.5830   | 0.1228 | 0.2651 |
| 0.1 | 0.0620   | 0.5890   | 0.1253 | 0.2646 |
| 0.2 | 0.0472   | 0.6580   | 0.1068 | 0.2317 |
| 0.3 | 0.0443   | 0.8100   | 0.0948 | 0.1686 |
| 0.4 | 0.0293   | 0.9300   | 0.0719 | 0.1326 |
| 0.5 | 0.0325   | 0.9350   | 0.0868 | 0.1388 |
| 0.6 | 0.0367   | 0.9100   | 0.0869 | 0.1550 |
| 0.7 | 0.0561   | 0.8570   | 0.1370 | 0.1948 |
| 0.8 | 0.0536   | 0.7890   | 0.1113 | 0.1997 |
| 0.9 | 0.0743   | 0.6510   | 0.1502 | 0.2347 |

### 2a.: `hc_dist = "complete"`

| rho | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|-----|----------|----------|--------|--------|
| 0.0 | 0.1178   | 0.9090   | 0.1470 | 0.1791 |
| 0.1 | 0.1225   | 0.9020   | 0.1564 | 0.1870 |
| 0.2 | 0.1175   | 0.9230   | 0.1368 | 0.1469 |
| 0.3 | 0.0947   | 0.9300   | 0.1327 | 0.1428 |
| 0.4 | 0.0944   | 0.9420   | 0.1311 | 0.1091 |
| 0.5 | 0.0933   | 0.9340   | 0.1337 | 0.1171 |
| 0.6 | 0.0959   | 0.9170   | 0.1314 | 0.1349 |
| 0.7 | 0.1027   | 0.8650   | 0.1453 | 0.1628 |
| 0.8 | 0.1046   | 0.7800   | 0.1514 | 0.1965 |
| 0.9 | 0.1078   | 0.6500   | 0.1635 | 0.2257 |

### 2a.: `hc_dist = "average"`

| rho | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|-----|----------|----------|--------|--------|
| 0.0 | 0.1166   | 0.9110   | 0.1422 | 0.1756 |
| 0.1 | 0.1211   | 0.9040   | 0.1527 | 0.1837 |
| 0.2 | 0.1188   | 0.9230   | 0.1397 | 0.1483 |
| 0.3 | 0.0908   | 0.9310   | 0.1280 | 0.1426 |
| 0.4 | 0.0905   | 0.9430   | 0.1274 | 0.1087 |
| 0.5 | 0.0854   | 0.9390   | 0.1307 | 0.1173 |
| 0.6 | 0.0913   | 0.9120   | 0.1327 | 0.1384 |
| 0.7 | 0.1007   | 0.8550   | 0.1485 | 0.1724 |
| 0.8 | 0.1051   | 0.7770   | 0.1545 | 0.1987 |
| 0.9 | 0.1078   | 0.6500   | 0.1635 | 0.2239 |

---

## 3. AR(1) Blocks + White Block: Q sweep

We consider $n=150$, $p=1000$, $M=5$, $s=5$, $rho=0.7$, $SNR=2.0$, $200$ MCs,
$Q \in \{5, 10, 15, 20, 25, 30, 35, 40, 45, 50\}$
$p_{\text{ar}} = M * Q$ varies from 25 to 250 and $p_{\text{white}} = p - p_{\text{ar}}$.

### 3a.: `hc_dist = "single"`

| Q  | p_ar | p_white | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|------|---------|----------|----------|--------|--------|
| 5  | 25   | 975     | 0.0561   | 0.8570   | 0.1370 | 0.1948 |
| 10 | 50   | 950     | 0.0345   | 0.8560   | 0.0833 | 0.1935 |
| 15 | 75   | 925     | 0.0261   | 0.8210   | 0.0751 | 0.2539 |
| 20 | 100  | 900     | 0.0348   | 0.7860   | 0.1124 | 0.3044 |
| 25 | 125  | 875     | 0.0331   | 0.7230   | 0.1008 | 0.3562 |
| 30 | 150  | 850     | 0.0196   | 0.7510   | 0.0753 | 0.3370 |
| 35 | 175  | 825     | 0.0334   | 0.6540   | 0.0864 | 0.3663 |
| 40 | 200  | 800     | 0.0383   | 0.7170   | 0.0976 | 0.3069 |
| 45 | 225  | 775     | 0.0466   | 0.7420   | 0.1176 | 0.3180 |
| 50 | 250  | 750     | 0.0603   | 0.7380   | 0.1445 | 0.2594 |

### 3a.: `hc_dist = "complete"`

| Q  | p_ar | p_white | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|------|---------|----------|----------|--------|--------|
| 5  | 25   | 975     | 0.1027   | 0.8650   | 0.1453 | 0.1628 |
| 10 | 50   | 950     | 0.1026   | 0.8770   | 0.1308 | 0.1674 |
| 15 | 75   | 925     | 0.0733   | 0.8830   | 0.1136 | 0.1463 |
| 20 | 100  | 900     | 0.1020   | 0.8690   | 0.1414 | 0.1833 |
| 25 | 125  | 875     | 0.0994   | 0.8550   | 0.1527 | 0.1804 |
| 30 | 150  | 850     | 0.0840   | 0.8970   | 0.1326 | 0.1432 |
| 35 | 175  | 825     | 0.1032   | 0.8480   | 0.1399 | 0.1807 |
| 40 | 200  | 800     | 0.1004   | 0.8630   | 0.1408 | 0.1833 |
| 45 | 225  | 775     | 0.0951   | 0.8560   | 0.1340 | 0.1715 |
| 50 | 250  | 750     | 0.1227   | 0.8460   | 0.1610 | 0.1951 |

### 3a.: `hc_dist = "average"`

| Q  | p_ar | p_white | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|------|---------|----------|----------|--------|--------|
| 5  | 25   | 975     | 0.1007   | 0.8550   | 0.1485 | 0.1724 |
| 10 | 50   | 950     | 0.0988   | 0.8710   | 0.1278 | 0.1676 |
| 15 | 75   | 925     | 0.0755   | 0.8690   | 0.1281 | 0.1661 |
| 20 | 100  | 900     | 0.0983   | 0.8660   | 0.1369 | 0.1861 |
| 25 | 125  | 875     | 0.0976   | 0.8540   | 0.1502 | 0.1834 |
| 30 | 150  | 850     | 0.0857   | 0.8910   | 0.1422 | 0.1537 |
| 35 | 175  | 825     | 0.1114   | 0.8390   | 0.1487 | 0.1867 |
| 40 | 200  | 800     | 0.1084   | 0.8540   | 0.1484 | 0.1888 |
| 45 | 225  | 775     | 0.1081   | 0.8460   | 0.1581 | 0.1812 |
| 50 | 250  | 750     | 0.1233   | 0.8460   | 0.1573 | 0.1920 |

---

## 4. AR(1) Block + White Block: M sweep

We consider $n=150$, $p=1000$, $Q=5$, $rho=0.7$, $SNR=2.0$, $200$ MCs,
$M \in \{1, 3, 5, 10, 15, 20, 30\}$.
$p_{\text{ar}} = M * 5$ varies from 5 to 150 and $p_{\text{white}} = p - p_{\text{ar}}$,
$s = M$.

### 4a.: `hc_dist = "single"`

| M  | s  | p_ar | p_white | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|----|------|---------|----------|----------|--------|--------|
| 1  | 1  | 5    | 995     | 0.0100   | 0.9000   | 0.0702 | 0.3008 |
| 3  | 3  | 15   | 985     | 0.0449   | 0.9017   | 0.1193 | 0.1854 |
| 5  | 5  | 25   | 975     | 0.0561   | 0.8570   | 0.1370 | 0.1948 |
| 10 | 10 | 50   | 950     | 0.0627   | 0.3525   | 0.1341 | 0.2064 |
| 15 | 15 | 75   | 925     | 0.0519   | 0.1013   | 0.1595 | 0.1201 |
| 20 | 20 | 100  | 900     | 0.0673   | 0.0380   | 0.2035 | 0.0600 |
| 30 | 30 | 150  | 850     | 0.0642   | 0.0093   | 0.2078 | 0.0225 |

### 4a.: `hc_dist = "complete"`

| M  | s  | p_ar | p_white | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|----|------|---------|----------|----------|--------|--------|
| 1  | 1  | 5    | 995     | 0.0317   | 0.9100   | 0.1361 | 0.2869 |
| 3  | 3  | 15   | 985     | 0.0703   | 0.9200   | 0.1351 | 0.1540 |
| 5  | 5  | 25   | 975     | 0.1027   | 0.8650   | 0.1453 | 0.1628 |
| 10 | 10 | 50   | 950     | 0.1143   | 0.3750   | 0.1617 | 0.2131 |
| 15 | 15 | 75   | 925     | 0.1054   | 0.1187   | 0.1973 | 0.1243 |
| 20 | 20 | 100  | 900     | 0.0988   | 0.0418   | 0.2367 | 0.0630 |
| 30 | 30 | 150  | 850     | 0.0885   | 0.0107   | 0.2453 | 0.0271 |

### 4a.: `hc_dist = "average"`

| M  | s  | p_ar | p_white | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|----|----|------|---------|----------|----------|--------|--------|
| 1  | 1  | 5    | 995     | 0.0267   | 0.9100   | 0.1175 | 0.2869 |
| 3  | 3  | 15   | 985     | 0.0624   | 0.9083   | 0.1230 | 0.1565 |
| 5  | 5  | 25   | 975     | 0.1007   | 0.8550   | 0.1485 | 0.1724 |
| 10 | 10 | 50   | 950     | 0.1045   | 0.3665   | 0.1440 | 0.2101 |
| 15 | 15 | 75   | 925     | 0.0979   | 0.1197   | 0.1921 | 0.1276 |
| 20 | 20 | 100  | 900     | 0.0938   | 0.0410   | 0.2279 | 0.0624 |
| 30 | 30 | 150  | 850     | 0.0786   | 0.0108   | 0.2240 | 0.0269 |

---
