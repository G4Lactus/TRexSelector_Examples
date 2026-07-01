# Memory-Mapped Matrix Utilities: Demonstration Suite

## Overview

This folder contains C++ examples for the **`MemoryMappedMatrix<Scalar>`** utility from `utils/memmap/`.

The utility provides a disk-backed matrix interface that can be used like an Eigen matrix through `Eigen::Map`, making it useful when matrices are too large to keep entirely in RAM. In the TRexSelector project, this is especially relevant for large design matrices and dummy matrices that exceed available RAM in solver and selector workflows.

The main goals of this folder are:

1. to show how memory-mapped matrices are created and accessed in C++,
2. to demonstrate out-of-core workflows for large matrices,
3. to illustrate safe file lifetime management via RAII,
4. to provide a C++ companion to the corresponding R example.

The demo in this folder mirrors `R/memory_mapping/demo_memory_mapping.R` and also showcases C++-specific features such as RAII, zero-copy access through `Eigen::Map`, and OpenMP-based parallel generation.

---

## What this folder covers

A memory-mapped matrix stores its entries in a file on disk and exposes them through a matrix-like interface:

$$
\boldsymbol{X}_{\mathrm{mmap}} \in \mathbb{R}^{n \times p}.
$$

Conceptually, the matrix behaves like a regular numeric matrix, but its storage is file-backed rather than purely in-memory.

This is useful for:

- **large matrices**, when RAM is limited,
- **persistent intermediate data**, when matrix contents should survive beyond a single object lifetime,
- **solver interoperability**, because the mapped data can be consumed through Eigen-compatible interfaces.

---

## Start here

If you are new to memory mapping in TRexSelector, begin with:

1. **`demo_memory_mapping/`**  
   A self-contained demo covering matrix creation, read/write access, serial and parallel out-of-core generation, and safe cleanup.

---

## Demo overview

| Folder                    | Purpose |
|---------------------------|---------|
| `demo_memory_mapping/`    | Demonstrates basic mmap usage, out-of-core writes, OpenMP-based parallel filling, and element-wise access checks |

---

## Building and running

From the C++ workspace root:

```bash
# navigate to the TRexSelector C++ workspace
cd TRexSelector_Simulations/cpp

# build the demo in debug mode, linking against the installed TRexSelector library
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"

# build the demo target
cmake --build build/debug --target demo_memory_mapping
```

Run the demo with:

```bash
./build/debug/bin/memory_mapping/demo_memory_mapping/demo_memory_mapping
```

---

## Folder contents

```txt
memory_mapping/
  ├── README.md
  ├── CMakeLists.txt
  └── demo_memory_mapping/
      ├── demo_memory_mapping.cpp
      └── README.md
```

---

## Notes for new users

- This demo prints results to the console and does not write simulation summary files.
- Temporary backing files are created in the system temporary directory and removed automatically.
- The demo is intended to explain the memory-mapping utility itself, not to benchmark maximum scale production performance.

---

**Last updated**: 2026-07-01
