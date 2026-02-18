#!/usr/bin/env Rscript
# ============================================================
# Manual test script for ScSpatialFeaturePlot
# Run from the repo root:
#   Rscript tests/run_tests.R
# ============================================================

pass <- 0L
fail <- 0L

check <- function(label, expr) {
  result <- tryCatch(
    {
      force(expr)
      TRUE
    },
    error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      FALSE
    }
  )
  if (isTRUE(result)) {
    cat(sprintf("  PASS  %s\n", label))
    pass <<- pass + 1L
  } else {
    cat(sprintf("  FAIL  %s\n", label))
    fail <<- fail + 1L
  }
}

cat("=== Loading package ===\n")
check("Package loads without error", {
  library(ScSpatialFeaturePlot)
})

# ----------------------------------------------------------
cat("\n=== GetVisiumData ===\n")
spat <- NULL
check("GetVisiumData returns a Seurat object", {
  spat <<- GetVisiumData()
  stopifnot(inherits(spat, "Seurat"))
})

if (is.null(spat)) {
  cat("\nCannot continue without example data. Stopping.\n")
  quit(status = 1)
}

genes <- rownames(spat)
gene1 <- genes[1]
gene2 <- if (length(genes) >= 2) genes[2] else gene1
gene3 <- if (length(genes) >= 3) genes[3] else gene1

cat(sprintf("\n  Using genes: %s, %s, %s\n", gene1, gene2, gene3))
cat(sprintf("  Cells: %d   Genes: %d\n", ncol(spat), nrow(spat)))

# ----------------------------------------------------------
cat("\n=== Input validation ===\n")

check("Rejects non-Seurat input", {
  tryCatch(
    { ScSpatialFeaturePlot(spat = "bad", features = "x"); stop("should have errored") },
    error = function(e) {
      if (!grepl("Seurat", conditionMessage(e))) stop(conditionMessage(e))
    }
  )
})

check("Rejects empty features vector", {
  tryCatch(
    { ScSpatialFeaturePlot(spat = spat, features = character(0)); stop("should have errored") },
    error = function(e) {
      if (!grepl("non-empty", conditionMessage(e))) stop(conditionMessage(e))
    }
  )
})

check("Rejects invalid flip value", {
  tryCatch(
    { ScSpatialFeaturePlot(spat = spat, features = gene1, flip = 99); stop("should have errored") },
    error = function(e) {
      if (!grepl("flip", conditionMessage(e))) stop(conditionMessage(e))
    }
  )
})

check("Rejects uq < lq", {
  tryCatch(
    { ScSpatialFeaturePlot(spat = spat, features = gene1, uq = 0.1, lq = 0.9); stop("should have errored") },
    error = function(e) {
      if (!grepl("uq", conditionMessage(e))) stop(conditionMessage(e))
    }
  )
})

# ----------------------------------------------------------
cat("\n=== Single feature plot ===\n")

check("Returns ggplot for one gene", {
  p <- ScSpatialFeaturePlot(spat, features = gene1)
  stopifnot(inherits(p, "ggplot"))
})

# ----------------------------------------------------------
cat("\n=== Multiple feature plot ===\n")

check("Returns ggplot for two genes", {
  p <- ScSpatialFeaturePlot(spat, features = c(gene1, gene2))
  stopifnot(inherits(p, "ggplot"))
})

# ----------------------------------------------------------
cat("\n=== combine = FALSE ===\n")

check("Returns named list of ggplots", {
  feats <- c(gene1, gene2, gene3)
  plots <- ScSpatialFeaturePlot(spat, features = feats, combine = FALSE)
  stopifnot(is.list(plots))
  stopifnot(length(plots) == length(feats))
  stopifnot(all(names(plots) == feats))
  for (p in plots) stopifnot(inherits(p, "ggplot"))
})

# ----------------------------------------------------------
cat("\n=== Log transform ===\n")

check("logs = TRUE produces ggplot", {
  p <- ScSpatialFeaturePlot(spat, features = gene1, logs = TRUE)
  stopifnot(inherits(p, "ggplot"))
})

# ----------------------------------------------------------
cat("\n=== Flip modes ===\n")

for (fv in 1:3) {
  check(sprintf("flip = %d produces ggplot", fv), {
    p <- ScSpatialFeaturePlot(spat, features = gene1, flip = fv)
    stopifnot(inherits(p, "ggplot"))
  })
}

# ----------------------------------------------------------
cat("\n=== Custom parameters ===\n")

check("Custom pt.size.factor works", {
  p <- ScSpatialFeaturePlot(spat, features = gene1, pt.size.factor = 3)
  stopifnot(inherits(p, "ggplot"))
})

check("Custom quantile clipping works (uq=0.95, lq=0.05)", {
  p <- ScSpatialFeaturePlot(spat, features = gene1, uq = 0.95, lq = 0.05)
  stopifnot(inherits(p, "ggplot"))
})

check("Custom legend label works", {
  p <- ScSpatialFeaturePlot(spat, features = gene1, LegLabel = "Expression")
  stopifnot(inherits(p, "ggplot"))
})

check("ncol/nrow layout works", {
  p <- ScSpatialFeaturePlot(spat, features = c(gene1, gene2), ncol = 1)
  stopifnot(inherits(p, "ggplot"))
})

# ----------------------------------------------------------
cat("\n=== No mutation of input object ===\n")

check("Seurat object metadata unchanged after plotting", {
  meta_before <- colnames(spat[[]])
  ScSpatialFeaturePlot(spat, features = gene1)
  meta_after <- colnames(spat[[]])
  stopifnot(identical(meta_before, meta_after))
})

# ----------------------------------------------------------
cat("\n=== Saving a plot to file ===\n")

check("Plot saves to PNG without error", {
  p <- ScSpatialFeaturePlot(spat, features = c(gene1, gene2))
  tmp <- tempfile(fileext = ".png")
  ggplot2::ggsave(tmp, p, width = 8, height = 4, dpi = 72)
  stopifnot(file.exists(tmp))
  cat(sprintf("    Saved to: %s (%s bytes)\n", tmp, file.size(tmp)))
})

# ----------------------------------------------------------
cat("\n=== Unit tests (testthat) ===\n")
if (requireNamespace("testthat", quietly = TRUE)) {
  test_results <- testthat::test_dir(
    system.file("tests", "testthat", package = "ScSpatialFeaturePlot",
                mustWork = FALSE),
    reporter = "summary"
  )
  # If tests dir isn't installed, run from source
  if (is.null(test_results) || length(test_results) == 0) {
    test_dir <- file.path(getwd(), "tests", "testthat")
    if (dir.exists(test_dir)) {
      cat("  Running testthat from source tree...\n")
      testthat::test_dir(test_dir, reporter = "summary")
    } else {
      cat("  testthat directory not found, skipping.\n")
    }
  }
} else {
  cat("  testthat not installed, skipping.\n")
}

# ----------------------------------------------------------
cat("\n============================\n")
cat(sprintf("Results: %d passed, %d failed\n", pass, fail))
cat("============================\n")

if (fail > 0) quit(status = 1)
