# ==============================================================================
# demo_memory_mapping.py
# ==============================================================================
#
# Demonstrates the memory-mapping utilities provided by the trex_selector
# package.  Memory-mapped matrices allow handling large datasets without
# loading the entire array into Python's active memory, seamlessly passing
# a pointer to the C++ backend.
#
# Mirrors R/memory_mapping/demo_memory_mapping.R.
#
# ==============================================================================

import os
import tempfile

import numpy as np
from trex_selector.utils import AccessMode, MemoryMappedMatrix, numpy_to_memmap

print("=" * 65)
print(" Memory Mapping Demo for TRexSelector")
print("=" * 65 + "\n")

# ==============================================================================
# 1. Create a sample in-memory matrix
# ==============================================================================

n_rows = 100
n_cols = 50
print(f"1. Creating a sample {n_rows} x {n_cols} numeric matrix.\n")

rng = np.random.default_rng(42)
X = rng.standard_normal((n_rows, n_cols))


# ==============================================================================
# 2. Define a temporary file path
# ==============================================================================

fd, bin_file = tempfile.mkstemp(suffix=".bin")
os.close(fd)   # release OS file descriptor; numpy_to_memmap will re-open the path
print(f"2. Temporary file for memory mapping: {bin_file}\n")


# ==============================================================================
# Main demo — sections 3-8 wrapped in try/finally for safe cleanup
# ==============================================================================

try:
    # ==========================================================================
    # 3. Convert to Memory-Mapped Matrix
    # ==========================================================================

    print("3. Converting the in-memory matrix to a MemoryMappedMatrix ...")
    X_mmap = numpy_to_memmap(bin_file, X)
    print("   Done.\n")


    # ==========================================================================
    # 4. Explore the object
    # ==========================================================================

    print("4. Exploring the MemoryMappedMatrix object:")
    print(f"   rows()  = {X_mmap.rows()}")
    print(f"   cols()  = {X_mmap.cols()}")
    print(f"   size()  = {X_mmap.size()}   (rows * cols = {X_mmap.rows() * X_mmap.cols()})\n")


    # ==========================================================================
    # 5. Extract a block (subsetting via the numpy view)
    # ==========================================================================

    print("5. Extracting a 3x3 block via to_numpy()[0:3, 0:3]:")
    block = X_mmap.to_numpy()[0:3, 0:3]
    print(np.round(block, 4))
    match = np.allclose(block, X[0:3, 0:3])
    print(f"   Matches original matrix block? {match}\n")


    # ==========================================================================
    # 6. Read back into memory
    # ==========================================================================

    print("6. Reading the full memory-mapped matrix back into a NumPy array:")
    # .to_numpy() returns a zero-copy view; .copy() gives an independent array.
    X_restored = X_mmap.to_numpy().copy()
    match_full = np.allclose(X, X_restored)
    print(f"   Matches original full matrix? {match_full}\n")


    # ==========================================================================
    # 7. Re-open the file in read-only mode
    # ==========================================================================

    print("7. Re-opening the existing binary file in ReadOnly mode:")
    X_ro = MemoryMappedMatrix(bin_file, n_rows, n_cols, AccessMode.ReadOnly)
    print(f"   rows() = {X_ro.rows()},  cols() = {X_ro.cols()}\n")


    # ==========================================================================
    # 8. Out-of-core data generation: write column-by-column to a new file
    # ==========================================================================

    print("8. Out-of-core generation: writing data column-by-column to a new file.\n")
    print("   For large datasets, streaming data column-by-column to disk is")
    print("   significantly more memory-efficient than building the entire array in")
    print("   RAM first. Element-wise access via mmap[i, j] is demonstrated in")
    print("   Section 9 below.\n")

    fd_ooc, bin_file_ooc = tempfile.mkstemp(suffix=".bin")
    os.close(fd_ooc)

    try:
        rng_ooc = np.random.default_rng(0)
        print("   Writing data column-by-column to stream it efficiently ...")
        with open(bin_file_ooc, "wb") as fh:
            for _ in range(n_cols):
                # Each column is written as a contiguous float64 block.
                # Column-major (Fortran) order matches C++ Eigen ColMajor storage.
                col = rng_ooc.standard_normal(n_rows).astype(np.float64)
                fh.write(col.tobytes())

        print("   Mapping the sequentially populated binary file ...")
        X_ooc = MemoryMappedMatrix(bin_file_ooc, n_rows, n_cols, AccessMode.ReadOnly)
        print(f"   rows() = {X_ooc.rows()},  cols() = {X_ooc.cols()}")
        print(f"   First element: {X_ooc[0, 0]:.6f}   (read via __getitem__ — to_numpy() requires ReadWrite mode)\n")
    finally:
        try:
            os.unlink(bin_file_ooc)
            print("   -> Out-of-core temporary file cleaned up.")
        except OSError:
            pass


    # ==========================================================================
    # 9. Element-wise access: __getitem__ and __setitem__
    # ==========================================================================

    print("9. Element-wise access via mmap[row, col] (0-based indices).\n")

    # --- 9a. Scalar read ---
    r, c = 2, 3
    val = X_mmap[r, c]
    assert isinstance(val, float), "Expected a plain float scalar"
    match_scalar = np.isclose(val, X_mmap.to_numpy()[r, c])
    print(f"   Read  X_mmap[{r}, {c}] = {val:.6f}   (matches numpy view: {match_scalar})")

    # --- 9b. Scalar write ---
    original_val = X_mmap[r, c]
    X_mmap[r, c] = 99.0
    readback = X_mmap[r, c]
    numpy_readback = X_mmap.to_numpy()[r, c]
    print(f"   Write X_mmap[{r}, {c}] = 99.0 -> read back = {readback:.1f}   "
          f"(numpy view agrees: {np.isclose(numpy_readback, 99.0)})")

    # Revert to the original value so the rest of the demo is unaffected
    X_mmap[r, c] = original_val
    print(f"   Reverted to original value: {X_mmap[r, c]:.6f}")

    # --- 9c. Bounds guard ---
    print(f"   Bounds check (row index {n_rows} on a {n_rows}-row matrix):")
    try:
        _ = X_mmap[n_rows, 0]   # out-of-bounds
    except IndexError as e:
        print(f"   Caught IndexError: {e}")

    # --- 9d. Mode guard ---
    print(f"   Write attempt on ReadOnly matrix:")
    try:
        X_ro[0, 0] = 1.0        # X_ro was opened ReadOnly in Section 7
    except RuntimeError as e:
        print(f"   Caught RuntimeError: {e}")

    # --- 9e. Block write via slice __setitem__ ---
    block_rows = slice(0, 3)
    block_cols = slice(0, 4)
    new_block = np.arange(12.0, dtype=np.float64).reshape(3, 4)
    original_block = X_mmap.to_numpy()[block_rows, block_cols].copy()
    X_mmap[block_rows, block_cols] = new_block
    readback_block = X_mmap.to_numpy()[block_rows, block_cols]
    block_match = np.allclose(readback_block, new_block)
    print(f"   Block write X_mmap[0:3, 0:4] = arange(12).reshape(3,4)  "
          f"-> round-trip match: {block_match}")

    # Revert block
    X_mmap[block_rows, block_cols] = original_block
    print(f"   Reverted block.")

    # --- 9f. Explicit write_block() method ---
    new_block2 = np.full((3, 4), fill_value=7.0, dtype=np.float64)
    original_block2 = X_mmap.to_numpy()[0:3, 0:4].copy()
    X_mmap.write_block(0, 3, 0, 4, new_block2)
    readback_block2 = X_mmap.to_numpy()[0:3, 0:4]
    block2_match = np.allclose(readback_block2, 7.0)
    print(f"   write_block(0,3, 0,4, fill=7.0)            "
          f"-> round-trip match: {block2_match}")

    # Revert
    X_mmap.write_block(0, 3, 0, 4, original_block2)
    print(f"   Reverted block.")

    print()

finally:
    try:
        os.unlink(bin_file)
        print("   -> Temporary file cleaned up.")
    except OSError:
        pass

print("\nDemo completed successfully!")
