# Demo: Memory-Mapped Matrix Usage

## Purpose

This demo illustrates how to use `MemoryMappedMatrix` from Python.

It mirrors the corresponding C++ example in
[cpp/memory_mapping/demo_memory_mapping/](../../../cpp/memory_mapping/demo_memory_mapping/)
(and the R example in `R/memory_mapping/demo_memory_mapping.R`), replacing the
C++-specific features (RAII cleanup, `Eigen::Map` interoperability, OpenMP
parallel generation) with their Python equivalents: `try/finally` cleanup,
zero-copy NumPy views via `to_numpy()`, and deterministic per-column seeding.

The demo is organized into four scenarios:

1. **Basics**: create a memory-mapped matrix, write data into it, inspect
   metadata, extract a block, verify a full roundtrip, and reopen the file in
   read-only mode.
2. **Out-of-core streaming (serial)**: generate a matrix column by column
   without allocating the full matrix in RAM.
3. **Per-column seeded fill**: fill columns through the zero-copy view with
   deterministic per-column random seeds. The OpenMP-parallel version of this
   loop is C++-exclusive; the seeding pattern that makes it deterministic is
   reproduced here serially.
4. **Element-wise access**: test scalar reads, scalar writes, compound
   assignment, bounds checking, read-only protection, and the Python
   extensions (slice `__setitem__` block writes and `write_block()`).

---

## What this demo shows

The demo focuses on several practical ideas:

- how a file-backed matrix can be used like a NumPy array,
- how to reduce RAM usage by writing data directly into a mapped file,
- how `try/finally` blocks ensure that temporary files are cleaned up safely,
- how deterministic column-wise random number generation supports reproducible
  parallel workflows.

Entry points used (from `trex_selector_neo.utils`):

- `MemoryMappedMatrix` — disk-backed matrix class (`rows()`, `cols()`,
  `size()`, `to_numpy()`, `write_block()`, `__getitem__`/`__setitem__`)
- `numpy_to_memmap(path, X)` — write a NumPy array to a binary file and
  return a `MemoryMappedMatrix` over it
- `AccessMode` — `ReadOnly` / `ReadWrite` open modes

---

## Running the demo

```bash
python demo_memory_mapping.py
```

The demo writes its output to the console only. No `simulation_results/`
folder is created for this example.

---

## What to look for

When reading the console output, check the following points:

- the created matrix has the expected dimensions,
- extracted blocks and full roundtrip checks match the original in-memory data,
- read-only reopening works as expected,
- serial out-of-core generation reproduces the expected column values,
- per-column seeds reproduce each column independently of fill order,
- invalid access triggers the intended exception handling
  (`IndexError` for out-of-bounds, `RuntimeError` for read-only writes),
- temporary `.bin` backing files are removed automatically after use.

---

## Technical notes

- The backing files are created in the system temporary directory.
- Cleanup is managed through `try/finally` blocks; the mapped objects are
  deleted (unmapped) before the backing file is removed.
- Binary files store columns as contiguous float64 blocks (column-major),
  matching the Eigen ColMajor storage used by the C++ backend.
- `to_numpy()` returns a zero-copy view that requires `ReadWrite` mode;
  read-only matrices are accessed element-wise via `mmap[i, j]`.

---

**Last updated**: 2026-07-21
