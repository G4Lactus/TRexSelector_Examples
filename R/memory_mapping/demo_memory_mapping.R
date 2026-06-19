# ==============================================================================
# demo_memory_mapping.R
# ==============================================================================
#
# This script demonstrates the memory-mapping functionalities for the
# TRexSelector package.
# Memory-mapped matrices allow handling large datasets without loading
# the entire matrix into R's active memory, seamlessly passing a pointer to the
# C++ backend.
#
# Note: File handling in this demo
# Memory-mapped matrices are written to a temporary binary in a temporary
# directory. The file is automatically cleaned up when the R session ends.
#
# ==============================================================================

library(TRexSelector)

cat("=================================================================\n")
cat(" Memory Mapping Demo for TRexSelector\n")
cat("=================================================================\n\n")

# 1. Create a sample in-memory matrix
n_rows <- 100
n_cols <- 50
cat(sprintf("1. Creating a sample %d x %d numeric matrix in R.\n", n_rows, n_cols))
set.seed(42)
mat <- matrix(rnorm(n_rows * n_cols), nrow = n_rows, ncol = n_cols)


# 2. Define a temporary file path
bin_file <- tempfile(pattern = "trex_data_", fileext = ".bin")
cat(sprintf("2. Temporary file for memory mapping: %s\n", bin_file))

# Ensure cleanup on exit
on.exit({
  if (file.exists(bin_file)) {
    unlink(bin_file)
    cat("   -> Temporary file cleaned up.\n")
  }
}, add = TRUE)


# 3. Convert to Memory-Mapped Matrix
cat("\n3. Converting the in-memory matrix to a MemoryMappedMatrix...\n")
mmap_mat <- convert_to_memory_mapped(mat, bin_file)


# 4. Explore the object using S3 methods
cat("\n4. Exploring the mmap_matrix object:\n")
print(mmap_mat)

cat("\n   Dimensions of the mmap_matrix:\n")
print(dim(mmap_mat))


# 5. Extract blocks (subsetting)
cat("\n5. Extracting a 3x3 block using the '[' operator:\n")
block <- mmap_mat[1:3, 1:3]
print(block)

# Verify against the original matrix
cat("   Matches original matrix block? ", isTRUE(all.equal(block, mat[1:3, 1:3])), "\n")


# 6. Read back into memory
cat("\n6. Reading the full memory-mapped matrix back into an R matrix:\n")
mat_restored <- as.matrix(mmap_mat)
cat("   Matches original full matrix? ", isTRUE(all.equal(mat, mat_restored)), "\n")


# 7. Re-open the file in read-only mode
cat("\n7. Re-opening the existing binary file in 'readonly' mode:\n")
mmap_readonly <- mmap_matrix(bin_file, rows = n_rows, cols = n_cols, mode = "readonly")
print(mmap_readonly)


# 8. Out-of-core data generation: Create file on disk and fill it incrementally
cat("\n8. Out-of-core generation: Creating file on disk and filling it incrementally:\n")
cat("   For large datasets, streaming data column-by-column to disk is significantly\n")
cat("   more memory-efficient than building the entire matrix in RAM first.\n")
cat("   Element-wise access via [ and [<- is demonstrated in Section 9 below.\n\n")

bin_file_ooc <- tempfile(pattern = "trex_data_ooc_", fileext = ".bin")

# Ensure cleanup for the second file
on.exit({
  if (file.exists(bin_file_ooc)) {
    unlink(bin_file_ooc)
    cat("   -> Out-of-core temporary file cleaned up.\n")
  }
}, add = TRUE)

cat("   Writing data column-by-column to stream it efficiently...\n")
con <- file(bin_file_ooc, "wb")
for (i in seq_len(n_cols)) {
  writeBin(as.double(rnorm(n_rows)), con)
}
close(con)

cat("   Mapping the sequentially populated binary file...\n")
mmap_ooc <- mmap_matrix(bin_file_ooc, rows = n_rows, cols = n_cols)
print(mmap_ooc)


# 9. Element-wise access: [ (read) and [<- (write)
cat(strrep("=", 65), "\n")
cat("9. Element-wise access via [ and [<- (1-based indices).\n\n")

# --- 9a. Scalar read returns a plain double, not a 1x1 matrix ---
r <- 2L
cc <- 3L
val <- mmap_mat[r, cc]
cat(sprintf("   Read  mmap_mat[%d, %d] = %.6f\n", r, cc, val))
cat(sprintf("   Is plain scalar (not matrix): %s\n", !is.matrix(val)))
cat(sprintf("   Matches original mat[%d, %d]? %s\n\n", r, cc,
            isTRUE(all.equal(val, mat[r, cc]))))

# --- 9b. Scalar write ---
orig_val <- mmap_mat[r, cc]
mmap_mat[r, cc] <- 99.0
readback <- mmap_mat[r, cc]
cat(sprintf("   Write mmap_mat[%d, %d] <- 99.0  -> read back = %.1f   (correct: %s)\n\n",
            r, cc, readback, isTRUE(all.equal(readback, 99.0))))

# Revert
mmap_mat[r, cc] <- orig_val
cat(sprintf("   Reverted to original value: %.6f\n\n", mmap_mat[r, cc]))

# --- 9c. Block write ---
block_rows <- 1:3
block_cols <- 1:4
new_block  <- matrix(seq_len(length(block_rows) * length(block_cols)), nrow = length(block_rows))
mmap_mat[block_rows, block_cols] <- new_block
readback_block <- mmap_mat[block_rows, block_cols]
cat(sprintf("   Block write mmap_mat[1:3, 1:4] <- matrix(1:12):\n"))
cat(sprintf("   Round-trip matches? %s\n\n", isTRUE(all.equal(readback_block, new_block))))

# Revert block to original
mmap_mat[block_rows, block_cols] <- mat[block_rows, block_cols]

# --- 9d. Mode guard ---
cat("   Write attempt on ReadOnly matrix:\n")
tryCatch(
  mmap_readonly[1, 1] <- 0,
  error = function(e) cat(sprintf("   Caught error: %s\n", conditionMessage(e)))
)


cat("\nDemo completed successfully!\n")
