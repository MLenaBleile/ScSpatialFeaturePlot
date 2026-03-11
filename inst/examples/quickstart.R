#!/usr/bin/env Rscript
# ============================================================
# ScSpatialFeaturePlot — Quick-Start Installation & Usage Demo
# ============================================================
#
# Download and run this script to verify the package works:
#
#   Rscript quickstart.R
#
# Prerequisites: R (>= 4.0) with 'remotes' or 'devtools' installed.
# ============================================================

# --- 1. Install the package (skip if already installed) ------
if (!requireNamespace("ScSpatialFeaturePlot", quietly = TRUE)) {
  cat("Installing ScSpatialFeaturePlot from GitHub...\n")
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("MLenaBleile/ScSpatialFeaturePlot")
}

# --- 2. Load libraries ---------------------------------------
library(ScSpatialFeaturePlot)
library(ggplot2)
cat("ScSpatialFeaturePlot loaded successfully.\n\n")

# --- 3. Load bundled example Visium data ---------------------
spat <- GetVisiumData()
cat(sprintf("Example data: %d genes x %d cells\n",
            nrow(spat), ncol(spat)))
cat(sprintf("Available genes: %s\n\n",
            paste(head(rownames(spat), 10), collapse = ", ")))

# --- 4. Single gene plot -------------------------------------
cat("Plotting single gene...\n")
p1 <- ScSpatialFeaturePlot(spat, features = "MLANA", pt.size.factor = 2)
ggsave("demo_single_gene.png", p1, width = 6, height = 5, dpi = 150)
cat("  Saved: demo_single_gene.png\n\n")

# --- 5. Multi-gene comparison (shared scale) -----------------
cat("Plotting multiple genes on shared scale...\n")
p2 <- ScSpatialFeaturePlot(spat,
                            features = c("MLANA", "MITF"),
                            pt.size.factor = 2)
ggsave("demo_multi_gene.png", p2, width = 10, height = 5, dpi = 150)
cat("  Saved: demo_multi_gene.png\n\n")

# --- 6. Separate scales (combine = FALSE) --------------------
cat("Plotting genes on separate scales...\n")
plots <- ScSpatialFeaturePlot(spat,
                               features = c("MLANA", "MITF", "SOX10", "PTPRC"),
                               combine = FALSE)
combined <- do.call(gridExtra::grid.arrange, c(plots, ncol = 2))
ggsave("demo_separate_scales.png", combined, width = 10, height = 10, dpi = 150)
cat("  Saved: demo_separate_scales.png\n\n")

# --- 7. Custom options: quantile clipping & flip -------------
cat("Plotting with quantile clipping and flip...\n")
p3 <- ScSpatialFeaturePlot(spat,
                            features = "VIM",
                            pt.size.factor = 2,
                            lq = 0.25, uq = 0.95,
                            flip = 1)
ggsave("demo_clipped_flipped.png", p3, width = 6, height = 5, dpi = 150)
cat("  Saved: demo_clipped_flipped.png\n\n")

# --- Done! ---------------------------------------------------
cat("All demos completed successfully!\n")
cat("Output files:\n")
for (f in list.files(pattern = "^demo_.*\\.png$")) {
  cat(sprintf("  %s  (%s KB)\n", f, round(file.size(f) / 1024)))
}
