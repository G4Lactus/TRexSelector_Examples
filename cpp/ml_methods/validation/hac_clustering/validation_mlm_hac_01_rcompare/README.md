# Validation: Single-Linkage HAC vs. R `hclust`

## Purpose

This validation performs a head-to-head numerical equivalence check between the C++ single-linkage hierarchical clustering implementation and R's `hclust(method="single")`.

Both implementations are run on the **same** dumped design matrices $X$, and the validation checks whether they produce the same dendrogram geometry and the same flat clusterings. This is a heavier workflow-level validation, not a lightweight unit test.

---

## Reference setup

The companion R workflow clusters the columns of each dataset using

$$
d(i,j) = 1 - |\mathrm{cor}(x_i, x_j)|
$$

followed by

```r
hclust(as.dist(d), method = "single")
```

and writes reference outputs to disk.

The C++ side reads the same matrices $X_n$ from CSV, runs the single-linkage HAC engine on the columns of $X$, and compares the resulting dendrogram against the R output.

---

## What is compared

For each dataset, the validation performs two comparisons.

### 1. Merge heights

Let

$$
h_1, h_2, \dots, h_{p-1}
$$

denote the sorted merge heights from the dendrogram. The program compares the sorted C++ merge distances against the sorted R merge heights and reports

$$
\max_k |h_k^{\mathrm{cpp}} - h_k^{\mathrm{R}}|.
$$

This comparison is robust to differences in internal node numbering or merge-row ordering, since the heights themselves are the invariant geometric quantities.

### 2. Flat partitions for every $K$

For every cluster count

$$
K = 2, 3, \dots, p-1,
$$

the C++ dendrogram is cut into $K$ clusters and compared to R's `cutree(hc, k = K)` output.

The agreement is measured with the Adjusted Rand Index (ARI),

$$
\mathrm{ARI} = 1
$$

meaning the two partitions are identical up to relabeling. The validation reports the minimum ARI across all tested $K$ as well as the number of cluster counts for which $\mathrm{ARI} < 1$.

---

## Pass criterion

A dataset passes if both of the following hold:

- merge heights agree to within
  $$
  10^{-9},
  $$
- every flat partition satisfies
  $$
  \mathrm{ARI} = 1.
  $$

The full validation passes only if this is true for **all** datasets in the reference dump.

---

## Running

```bash
./build/debug/bin/ml_methods/validation/hac_clustering/validation_mlm_hac_01_rcompare/validation_mlm_hac_01_rcompare
```

By default, the program reads its reference CSV files from the hard-coded dump directory in the source, but a different location can be supplied via the command-line argument:

```bash
./build/debug/bin/ml_methods/validation/hac_clustering/validation_mlm_hac_01_rcompare/validation_mlm_hac_01_rcompare --dir /path/to/rdump_tlars
```

---

## Expected inputs

The validation expects R-generated reference files including:

- `meta.csv`,
- `Xn_<i>.csv`,
- `r_clust_height_<i>.csv`,
- `r_clust_labels_<i>.csv`.

Here, `r_clust_height_<i>.csv` stores the sorted merge heights, and `r_clust_labels_<i>.csv` stores the R cluster labels for all cut levels.

---

## What to look for

When reading the console output, focus on these quantities:

- `max|Δheight|`, the largest absolute difference in merge heights,
- `minARI(K)`, the worst ARI over all tested cluster counts,
- `#K ARI<1`, the number of cuts where the partitions disagree,
- the final dataset-level and global `PASS` or `FAIL` verdicts.

A perfect validation run should show height differences below tolerance, `minARI(K) = 1`, and zero partition mismatches for every dataset.

---

## Notes

- The validation uses one OpenMP thread in order to keep the comparison deterministic and easier to interpret.
- The comparison is performed on the columns of $X$, using the correlation-based distance $1 - |\mathrm{cor}|$.
- The current source describes the factor-model datasets as having distinct distances, which makes the single-linkage dendrogram effectively unique in this setting.

---

**Last updated**: 2026-07-01
