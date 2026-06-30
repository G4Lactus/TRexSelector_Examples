# r_dummy_variance_probe.R
# Quantify R trex()'s OWN dummy-RNG variance on the identical dumped X.
# Reuses rdump/X_*.csv + r_lambda2.csv + truth.csv (no DGP rerun).
# For several internal seeds, recompute R's mean FDR -> compare to the C++ band.

suppressMessages(library(TRexSelector))

get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}
rdump <- file.path(get_script_dir(), "rdump")

num_MC     <- 100
target_fdr <- 0.10
M          <- 3

# ---- load truth + lambda2 ----
tr  <- read.csv(file.path(rdump, "truth.csv"),     stringsAsFactors = FALSE)
lam <- read.csv(file.path(rdump, "r_lambda2.csv"), stringsAsFactors = FALSE)
truth <- lapply(tr$indices, function(s) as.integer(strsplit(s, "-")[[1]]) + 1L)  # ->1-based
names(truth) <- tr$mc
lam2 <- setNames(lam$lambda2, lam$mc)

run_once <- function(seed_off) {
  fdrs <- numeric(num_MC)
  for (mc in 0:(num_MC - 1)) {
    X  <- as.matrix(read.csv(file.path(rdump, sprintf("X_%d.csv", mc)), header = FALSE))
    sv <- svd(X)
    y1 <- (X %*% sv$v[, 1:M])[, 1]
    set.seed(20240601 + seed_off * 100003 + mc)   # vary R's internal dummy RNG
    rr <- trex(X = X, y = y1, tFDR = target_fdr, method = "trex+GVS",
               GVS_type = "EN", lambda_2_lars = lam2[[as.character(mc)]],
               verbose = FALSE)
    sel <- which(rr$selected_var > .Machine$double.eps)
    k   <- length(sel)
    tp  <- length(intersect(truth[[as.character(mc)]], sel))
    fdrs[mc + 1] <- if (k == 0) 0 else (k - tp) / k
  }
  mean(fdrs)
}

cat("R trex() dummy-RNG variance over 8 internal seeds (identical X + R lambda2):\n")
vals <- numeric(8)
for (s in 0:7) {
  vals[s + 1] <- run_once(s)
  cat(sprintf("  seed_off %d : mean FDR = %.4f\n", s, vals[s + 1]))
}
cat(sprintf("\n  R band: min=%.4f  max=%.4f  mean=%.4f\n", min(vals), max(vals), mean(vals)))
cat("  C++ band (seeds 1000..8000): 0.064 .. 0.076 (mean ~0.071)\n")
