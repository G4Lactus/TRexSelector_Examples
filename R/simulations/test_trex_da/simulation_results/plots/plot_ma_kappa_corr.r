# ==============================================================================
# plot_ma_kappa_corr.R
#
# Output:
#   corr_ma_kappa_panel.png  -- static PNG for markdown (via kaleido)
#   corr_ma_kappa_panel.html -- interactive Plotly (open in browser)
#
# One-time setup:
#   install.packages(c("plotly", "htmlwidgets"))
#   plotly::install_kaleido()
#
# DGP mirrors dgp_nn_snr.R exactly.
# ==============================================================================

library(plotly)
library(htmlwidgets)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
N      <- 300L
P      <- 1000L
RHO    <- 0.7
SEED   <- 2026L
P_VIS  <- 200L
KAPPAS <- c(1L, 3L, 5L, 10L)

W_PX   <- 1100L   # figure width  (pixels)
H_PX   <- 950L    # figure height (pixels)
SCALE  <- 2L      # retina / HiDPI multiplier for PNG

# ------------------------------------------------------------------------------
# MA(kappa) DGP — identical to dgp_nn_snr()
# ------------------------------------------------------------------------------
gen_ma_kappa <- function(n, p, kappa, rho, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  l_seq <- 0L:kappa
  theta  <- rho^l_seq
  theta  <- theta / sqrt(sum(theta^2))
  eta    <- matrix(rnorm(n * (p + kappa)), nrow = n, ncol = p + kappa)
  X      <- matrix(0.0, nrow = n, ncol = p)
  for (j in seq_len(p))
    X[, j] <- drop(eta[, j:(j + kappa)] %*% theta)
  X
}

# RdBu colorscale (blue = -1, white = 0, red = +1)
RDBU <- list(
  list(0.0, "#053061"), list(0.1, "#2166ac"), list(0.2, "#4393c3"),
  list(0.3, "#92c5de"), list(0.4, "#d1e5f0"), list(0.5, "#f7f7f7"),
  list(0.6, "#fddbc7"), list(0.7, "#f4a582"), list(0.8, "#d6604d"),
  list(0.9, "#b2182b"), list(1.0, "#67001f")
)

# ------------------------------------------------------------------------------
# Build one heatmap per kappa
# ------------------------------------------------------------------------------
make_heatmap <- function(kappa, show_colorbar = FALSE) {
  cat(sprintf("  Generating MA(%d) ...\n", kappa))
  X <- gen_ma_kappa(N, P, kappa, RHO, SEED)
  C <- cor(X[, seq_len(P_VIS)])

  plot_ly(
    x            = seq_len(P_VIS),
    y            = seq_len(P_VIS),
    z            = C,
    type         = "heatmap",
    zmin         = -1,
    zmax         =  1,
    colorscale   = RDBU,
    reversescale = FALSE,
    showscale    = show_colorbar,
    colorbar     = if (show_colorbar) list(
      title     = list(text = "Corr", side = "right"),
      thickness = 14,
      len       = 0.45,
      x         = 1.01,
      y         = 0.77
    ) else NULL
  ) |>
    layout(
      title = list(
        text    = sprintf("<b>\u03ba = %d</b>", kappa),
        font    = list(size = 16),
        x       = 0.5,
        xanchor = "center"
      ),
      xaxis = list(
        title    = "Predictor j",
        tickfont = list(size = 9)
      ),
      yaxis = list(
        title     = "Predictor k",
        tickfont  = list(size = 9),
        autorange = "reversed"
      )
    )
}

# ------------------------------------------------------------------------------
# Build 2x2 panel
# ------------------------------------------------------------------------------
cat("Building plots ...\n")
plots <- list(
  make_heatmap( 1L, show_colorbar = FALSE),
  make_heatmap( 3L, show_colorbar = FALSE),
  make_heatmap( 5L, show_colorbar = FALSE),
  make_heatmap(10L, show_colorbar = TRUE)
)

fig <- subplot(
  plots[[1]], plots[[2]], plots[[3]], plots[[4]],
  nrows  = 2,
  shareX = FALSE,
  shareY = FALSE,
  titleX = TRUE,
  titleY = TRUE,
  margin = 0.09
) |>
  layout(
    title = list(
      text = paste0(
        "<b>MA(\u03ba) Banded Correlation Structure</b>",
        "<br><sup>\u03ba \u2208 {1, 3, 5, 10}",
        " \u2022 \u03c1 = 0.7",
        " \u2022 n = 300, p = 1000",
        " \u2022 first 200 predictors shown</sup>"
      ),
      font    = list(size = 17),
      x       = 0.5,
      xanchor = "center"
    ),
    margin = list(t = 100, b = 60, l = 70, r = 20)
  )

# ------------------------------------------------------------------------------
# 1) Interactive HTML — open in any browser
# ------------------------------------------------------------------------------
html_file <- "corr_ma_kappa_panel.html"
saveWidget(fig, html_file, selfcontained = TRUE,
           title = "MA(kappa) Correlation Matrices")
cat(sprintf("Interactive HTML: %s\n", html_file))

# ------------------------------------------------------------------------------
# 2) Static PNG — width/height passed to save_image(), NOT to layout()
#    Requires kaleido: plotly::install_kaleido()
# ------------------------------------------------------------------------------
png_file <- "corr_ma_kappa_panel.png"
save_image(fig, png_file, width = W_PX, height = H_PX, scale = SCALE)
cat(sprintf("Static PNG:       %s  (%dx%d px, scale=%d)\n",
            png_file, W_PX * SCALE, H_PX * SCALE, SCALE))

cat("\nDone.\n")
cat("  Embed in Markdown with:  ![](corr_ma_kappa_panel.png)\n")
