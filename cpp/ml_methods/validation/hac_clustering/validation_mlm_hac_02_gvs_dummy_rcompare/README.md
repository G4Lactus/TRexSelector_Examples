# Validation: GVS Dummy Matrix Cluster/Covariance Consistency

## Purpose

This validation checks that the C++ GVS dummy-generation setup is based on the **same dendrogram-cut clusters** and the **same per-cluster generating covariance** as R's `TRexSelector:::add_dummies_GVS()`.

The goal is not to compare the realized dummy matrices entry by entry, since those are random draws. Instead, the program compares the two deterministic ingredients that fully define how the dummies are generated from the clustering result.

---

## Reference setup

On the R side, the grouped dummy-generation workflow is:

1. compute the correlation distance
   $$
   d(i,j) = 1 - |\mathrm{cor}(x_i, x_j)|,
   $$
2. run

   ```r
   hclust(d, method = "single")
   ```

3. cut the dendrogram at
   $$
   h = 1 - \mathrm{corr\_max},
   $$
4. generate dummy blocks from a multivariate normal distribution using the within-cluster covariance matrix.

The C++ side follows the same conceptual recipe: single-linkage clustering under the same correlation distance, cutting the dendrogram at the same height, and building per-cluster covariance blocks from the clustered variables.

---

## What is compared

The validation compares the following two deterministic objects for every dataset and every value of `corr_max`.

### 1. Cluster partition at the cut height

The dendrogram is cut at

$$
h = 1 - \mathrm{corr\_max}.
$$

The C++ labels are obtained from `cut_tree_by_height(...)`, and the R labels come from `cutree(..., h = 1 - corr_max)`.

Agreement between the two partitions is measured with the Adjusted Rand Index (ARI). A value of

$$
\mathrm{ARI} = 1
$$

means that the partitions are identical up to relabeling.

### 2. Generating covariance for dummy blocks

For each cluster, the dummy variables are generated from a covariance matrix derived from the corresponding submatrix of $X$. The validation encodes this as a block-diagonal matrix

$$
\Sigma_{\mathrm{gen}},
$$

where

$$
\Sigma_{\mathrm{gen}}[i,j]
=
\begin{cases}
\mathrm{cov}(X)_{ij}, & \text{if variables } i \text{ and } j \text{ belong to the same cluster}, \\
0, & \text{otherwise}.
\end{cases}
$$

This is the mathematical covariance structure used to generate the grouped dummy blocks.

The program compares the C++ and R versions of this matrix through

$$
\max_{i,j} \left| \Sigma_{\mathrm{gen},ij}^{\mathrm{cpp}} - \Sigma_{\mathrm{gen},ij}^{\mathrm{R}} \right|.
$$


---

## Why the dummy draws themselves are not compared

The actual dummy matrices are random multivariate normal draws, so two correct implementations do not need to generate identical numeric samples. What matters is that they are sampled from the **same cluster structure** and the **same generating covariance**.

This validation therefore checks the deterministic setup behind the draws, which is the correct object to compare across languages.

---

## Pass criterion

For every pair `(dataset, corr_max)`, the validation requires:

- partition agreement with
  $$
  \mathrm{ARI} = 1,
  $$
- covariance agreement with
  $$
  \max |\Delta \Sigma_{\mathrm{gen}}| \le 10^{-9}.
  $$

The overall validation passes only if both conditions hold for every tested dataset and every value in the `corr_max` grid.

The source also notes that the C++ implementation may add a tiny ridge term of size $10^{-10}$ for Cholesky stability, but this is intentionally excluded from the compared $\Sigma_{\mathrm{gen}}$ object because the validation targets the mathematical generating covariance itself.

---

## Running

```bash
./build/debug/bin/ml_methods/validation/hac_clustering/validation_mlm_hac_02_gvs_dummy_rcompare/validation_mlm_hac_02_gvs_dummy_rcompare
```

As with `_01_`, the program reads reference CSV dumps from the default directory embedded in the source, but an alternative directory can be supplied via:

```bash
./build/debug/bin/ml_methods/validation/hac_clustering/validation_mlm_hac_02_gvs_dummy_rcompare/validation_mlm_hac_02_gvs_dummy_rcompare --dir /path/to/rdump_tlars
```

---

## Expected inputs

The validation expects R-generated reference files including:

- `meta.csv`,
- `gvs_corrmax.csv`,
- `Xn_<i>.csv`,
- `r_gvs_labels_<i>_<tag>.csv`,
- `r_gvs_sigma_<i>_<tag>.csv`.

Here, `<tag>` is the zero-padded encoding of the corresponding `corr_max` value, for example `050` for `0.50`.

---

## What to look for

When reading the console output, focus on these quantities:

- `corr_max`, the correlation threshold parameter,
- `h = 1 - corr_max`, the dendrogram cut height,
- `M`, the number of clusters obtained at that cut,
- `ARI(part)`, the partition agreement between C++ and R,
- `max|ΔΣ_gen|`, the maximum absolute covariance difference,
- the row-wise and final `PASS` or `FAIL` verdicts.

A perfect validation run should show `ARI(part) = 1` and covariance differences below tolerance for every dataset and every cut threshold.

---

## Notes

- The program runs with one OpenMP thread to keep the validation deterministic and easier to interpret.
- The compared partition is produced by cutting the single-linkage dendrogram at the inclusive threshold $h = 1 - \mathrm{corr\_max}$ on both the R and C++ sides.
- This validation is specifically about faithfulness of the **dummy-generation stage** in GVS, not about later selection outputs.

---

**Last updated**: 2026-07-01
