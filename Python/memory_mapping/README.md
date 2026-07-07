# Memory Mapping Demo (Python)

## Purpose

Demonstrate the memory-mapping utilities provided by the `trex_selector_neo`
package. Memory-mapped matrices allow handling large datasets without loading
the entire array into Python's active memory — the OS pages data in on demand,
and the same disk-backed buffer is passed as a pointer to the C++ backend.

Entry points used (from `trex_selector_neo.utils`):

- `MemoryMappedMatrix` — disk-backed matrix class (`rows()`, `cols()`,
  `size()`, `to_numpy()`, `write_block()`, `__getitem__`/`__setitem__`)
- `numpy_to_memmap(path, X)` — write a NumPy array to a binary file and
  return a `MemoryMappedMatrix` over it
- `AccessMode` — `ReadOnly` / `ReadWrite` open modes

---

## What the Demo Walks Through

[demo_memory_mapping.py](demo_memory_mapping.py) is a single script organized
into numbered sections:

1. **Create a sample in-memory matrix** — 100 x 50 standard-normal matrix.
2. **Define a temporary file path** — via `tempfile.mkstemp()`.
3. **Convert to Memory-Mapped Matrix** — `numpy_to_memmap(bin_file, X)`.
4. **Explore the object** — `rows()`, `cols()`, `size()`.
5. **Extract a block** — subsetting via the zero-copy `to_numpy()` view.
6. **Read back into memory** — `to_numpy().copy()` for an independent array,
   verified against the original.
7. **Re-open the file in read-only mode** —
   `MemoryMappedMatrix(path, n_rows, n_cols, AccessMode.ReadOnly)`.
8. **Out-of-core data generation** — stream data column-by-column to a new
   binary file (column-major float64 blocks, matching Eigen ColMajor storage),
   then map it read-only.
9. **Element-wise access** — scalar read/write via `mmap[i, j]`, bounds
   guard (`IndexError`), read-only write guard (`RuntimeError`), block write
   via slice `__setitem__`, and the explicit `write_block()` method.

Sections 3–9 are wrapped in `try/finally` blocks so temporary files are always
cleaned up.

---

## Running the Demo

```bash
python demo_memory_mapping.py
```

Output is printed to the console; all temporary binary files are removed
before the script exits. No `simulation_results/` folder is produced.

---

## Counterparts

- C++: [cpp/memory_mapping/demo_memory_mapping](../../cpp/memory_mapping/)
- R: `R/memory_mapping/demo_memory_mapping.R` (this script is a direct mirror)

---

**Last updated**: 2026-07-06
