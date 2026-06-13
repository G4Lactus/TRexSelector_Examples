# 2M AR(1) Blocks + White Noise

**Total:** $p = 2 M Q + p_\text{white}$.

### Covariance Matrix

$$
\boldsymbol{\Sigma} =
\begin{bmatrix}
\boldsymbol{\Sigma}_1(\rho) & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} \\
\mathbf{0} & \boldsymbol{\Sigma}_2(\rho) & \cdots & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} \\
\vdots & \vdots & \ddots & \vdots & \vdots &  & \vdots \\
\mathbf{0} & \mathbf{0} & \cdots & \boldsymbol{\Sigma}_M(\rho) & \mathbf{0} & \cdots & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} & \boldsymbol{\Sigma}_{M+1}(\rho) & \cdots & \mathbf{0} \\
\vdots & \vdots &  & \vdots & \vdots & \ddots & \vdots \\
\mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} & \boldsymbol{\Sigma}_{2M}(\rho) & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{I}_{p_\text{white}}
\end{bmatrix} \in \mathbb{R}^{p \times p}.
$$

### Support

$S = \{s_1,\dots,s_M\}$, $s_m \in \{(m-1)Q+1,\dots,mQ\}$, $\beta_{s_m} = \beta^* > 0$; blocks $M+1\dots2M$ + white = all nulls.

### Row Generation

$$
\boldsymbol{X}_i = \begin{pmatrix} 
\boldsymbol{X}_i^{(1)} \\ 
\vdots \\ 
\boldsymbol{X}_i^{(2M)} \\ 
\boldsymbol{X}_i^{(2M+1)} 
\end{pmatrix}, \quad \boldsymbol{X}_i^{(2M+1)} \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_{p_\text{white}}).
$$

### t-Extension ($\nu=3$)

$$
\boldsymbol{X}_i^{(m)} = \frac{\boldsymbol{Z}_i^{(m)}}{\sqrt{U_i^{(m)} / \nu}}, \quad m=1,\dots,2M+1,
$$

$\boldsymbol{Z}_i^{(m)} \sim \mathcal{N}(\mathbf{0}, \boldsymbol{\Sigma}_m(\rho))$ ($m=1\dots2M$), independent $U_i^{(m)}$.

**Cov:** $\text{Cov}(\boldsymbol{X}_i) = \frac{\nu}{\nu-2} \boldsymbol{\Sigma}$.

Penultimate diagonal: $\boldsymbol{\Sigma}_{2M}(\rho)$.
