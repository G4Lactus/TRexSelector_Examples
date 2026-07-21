# Memory-Mapped Matrix Utilities: Demonstration Suite

## Overview

This folder contains Python examples for the **`MemoryMappedMatrix`** utility
from `trex_selector_neo.utils`.

The utility provides a disk-backed matrix interface exposed to Python through
zero-copy NumPy views, making it useful when matrices are too large to keep
entirely in RAM. In the TRexSelector project, this is especially relevant for
large design matrices and dummy matrices that exceed available RAM in solver
and selector workflows — the same disk-backed buffer is passed as a pointer to
the C++ backend, with no copies on the way.

The main goals of this folder are:

1. to show how memory-mapped matrices are created and accessed from Python,
2. to demonstrate out-of-core workflows for large matrices,
3. to illustrate safe file lifetime management with `try/finally`.

The C++ demo suite additionally showcases C++-specific features such as RAII
and OpenMP-based parallel generation; see
[cpp/memory_mapping/](../../cpp/memory_mapping/) and `R/memory_mapping/` for
the corresponding walkthroughs.

---

## What this folder covers

A memory-mapped matrix stores its entries in a file on disk and exposes them
through a matrix-like interface:

$$
\boldsymbol{X}_{\mathrm{mmap}} \in \mathbb{R}^{n \times p}.
$$

Conceptually, the matrix behaves like a regular numeric matrix, but its
storage is file-backed rather than purely in-memory. The Python binding maps
`double` (float64) matrices.

This is useful for:

- **large matrices**, when RAM is limited,
- **persistent intermediate data**, when matrix contents should survive beyond
  a single object lifetime,
- **solver interoperability**, because the mapped data is consumed zero-copy
  by the C++ solvers and selectors.

---

## Start here

If you are new to memory mapping in TRexSelector, begin with:

1. **`demo_memory_mapping/`** — a self-contained demo covering matrix
   creation, read/write access, serial out-of-core generation, deterministic
   per-column seeding, and safe cleanup.

---

## Demo overview

| Folder | Purpose |
| ------ | ------- |
| `demo_memory_mapping/` | Demonstrates basic mmap usage, out-of-core writes, deterministic per-column seeded filling, and element-wise access checks |

---

## Running

```bash
python demo_memory_mapping/demo_memory_mapping.py
```

---

## Folder contents

```txt
memory_mapping/
  ├── README.md
  └── demo_memory_mapping/
      ├── demo_memory_mapping.py
      └── README.md
```

---

## Notes for new users

- This demo prints results to the console and does not write simulation
  summary files.
- Temporary backing files are created in the system temporary directory and
  removed automatically.
- The demo is intended to explain the memory-mapping utility itself, not to
  benchmark maximum scale production performance.

---

**Last updated**: 2026-07-21
