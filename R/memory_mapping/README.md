# Memory Mapping — R Demo

## Purpose

Demonstrates the memory-mapping functionality of the **TRexSelectorNeo**
package. Memory-mapped matrices allow handling large datasets without loading
the entire matrix into R's active memory: the R object seamlessly passes a
pointer to the disk-backed data through to the C++ backend.

---

## Demo

| File | Description |
|---|---|
| [demo_memory_mapping.R](demo_memory_mapping.R) | Create, access, fill, and clean up memory-mapped matrices via `convert_to_memory_mapped()` / `mmap_matrix()` |

The script walks through nine numbered sections:

1. Create a sample 100 x 50 in-memory matrix.
2. Define a temporary binary file path (auto-cleaned via `on.exit()`).
3. Convert the matrix to a `MemoryMappedMatrix` with `convert_to_memory_mapped()`.
4. Explore the `mmap_matrix` object via its S3 methods (`print`, `dim`).
5. Extract a 3 x 3 block with the `[` operator and verify against the original.
6. Read the full matrix back into an ordinary R matrix and verify equality.
7. Re-open the existing binary file in `"readonly"` mode with `mmap_matrix()`.
8. Out-of-core generation: create the file on disk and fill it column-by-column
   (streaming to disk instead of building the matrix in RAM first).
9. Element-wise access via `[` (read) and `[<-` (write), 1-based indices:
   scalar read (returns a plain double), scalar write, block write with
   round-trip check, and the write-guard error on a read-only matrix.

Memory-mapped files in this demo are written to a temporary directory and
removed when the R session ends.

---

## Running

```bash
Rscript R/memory_mapping/demo_memory_mapping.R
```

Output goes to the console only; no `simulation_results/` folder is produced.

---

## Counterparts

- C++: [cpp/memory_mapping/](../../cpp/memory_mapping/) — `demo_memory_mapping`
  (POSIX `mmap` + `Eigen::Map`, including OpenMP parallel filling).
- The same R mmap interface is exercised downstream by
  [../ml_methods/hac_clustering/demo_agg_hac_02_mmap.R](../ml_methods/hac_clustering/demo_agg_hac_02_mmap.R),
  [../ml_methods/standardization/demo_standardization_02_mmap.R](../ml_methods/standardization/demo_standardization_02_mmap.R),
  and the T-Rex demos 05–07 in
  [../trex_selector_methods/trex/](../trex_selector_methods/trex/README.md).

---

**Last updated**: 2026-07-06
