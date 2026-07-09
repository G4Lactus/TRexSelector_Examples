# Demo: Memory-Mapped Matrix Usage

## Purpose

This demo illustrates how to use `MemoryMappedMatrix<Scalar>` in C++.

It mirrors the corresponding R example in `R/memory_mapping/demo_memory_mapping.R`, but also includes C++-specific
features such as RAII-based cleanup, direct `Eigen::Map` interoperability, OpenMP-based parallel generation, and
element-wise access checks.

The demo is organized into four scenarios:

1. **Basics**: create a memory-mapped matrix, write data into it, inspect metadata, extract a block, verify a full
   roundtrip, and reopen the file in read-only mode.
2. **Out-of-core streaming (serial)**: generate a matrix column by column without allocating the full matrix in RAM.
3. **Out-of-core generation (parallel)**: fill columns in parallel with OpenMP using deterministic per-column random
   seeds.
4. **Element-wise access**: test scalar reads, scalar writes, compound assignment, bounds checking, and read-only
   protection.

---

## What this demo shows

The demo focuses on several practical ideas:

- how a file-backed matrix can be used like an Eigen object,
- how to reduce RAM usage by writing data directly into a mapped file,
- how C++ RAII helps ensure that temporary files are cleaned up safely,
- how deterministic column-wise random number generation supports reproducible parallel workflows.

---

## Running the demo

```bash
./build/debug/bin/memory_mapping/demo_memory_mapping/demo_memory_mapping
```

The demo writes its output to the console only. No `simulation_results/` folder is created for this example.

---

## What to look for

When reading the console output, check the following points:

- the created matrix has the expected dimensions,
- extracted blocks and full roundtrip checks match the original in-memory data,
- read-only reopening works as expected,
- serial out-of-core generation reproduces the expected column values,
- parallel generation remains deterministic through per-column seeds,
- invalid access triggers the intended exception handling,
- temporary `.bin` backing files are removed automatically after use.

---

## Technical notes

- The backing files are created in the system temporary directory.
- File cleanup is managed through an RAII guard so that unmapping happens before deletion.
- The parallel example uses OpenMP and assigns independent random seeds to columns, which avoids
  thread-scheduling-dependent randomness.

---

**Last updated**: 2026-07-01
