# ---- Input validation tests ----

test_that("ScSpatialFeaturePlot rejects non-Seurat input", {
  expect_error(
    ScSpatialFeaturePlot(spat = "not_a_seurat", features = "gene1"),
    "'spat' must be a Seurat object"
  )
  expect_error(
    ScSpatialFeaturePlot(spat = data.frame(), features = "gene1"),
    "'spat' must be a Seurat object"
  )
})

test_that("ScSpatialFeaturePlot rejects invalid features argument", {
  skip_if_not_installed("Seurat")
  # Create a minimal Seurat object
  counts <- matrix(rpois(200, 5), nrow = 10, ncol = 20)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell", 1:20)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  expect_error(
    ScSpatialFeaturePlot(spat = obj, features = 123),
    "'features' must be a non-empty character vector"
  )
  expect_error(
    ScSpatialFeaturePlot(spat = obj, features = character(0)),
    "'features' must be a non-empty character vector"
  )
})

test_that("ScSpatialFeaturePlot rejects invalid flip values", {
  skip_if_not_installed("Seurat")
  counts <- matrix(rpois(200, 5), nrow = 10, ncol = 20)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell", 1:20)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  expect_error(
    ScSpatialFeaturePlot(spat = obj, features = "Gene1", flip = 5),
    "'flip' must be 1, 2, or 3"
  )
  expect_error(
    ScSpatialFeaturePlot(spat = obj, features = "Gene1", flip = 0),
    "'flip' must be 1, 2, or 3"
  )
})

test_that("ScSpatialFeaturePlot rejects invalid quantile arguments", {
  skip_if_not_installed("Seurat")
  counts <- matrix(rpois(200, 5), nrow = 10, ncol = 20)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell", 1:20)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  expect_error(
    ScSpatialFeaturePlot(spat = obj, features = "Gene1", uq = 0.1, lq = 0.9),
    "'uq' and 'lq' must be numeric with uq >= lq"
  )
  expect_error(
    ScSpatialFeaturePlot(spat = obj, features = "Gene1", uq = "high"),
    "'uq' and 'lq' must be numeric"
  )
})

# ---- Functional tests (require spatial data) ----

test_that("ScSpatialFeaturePlot returns ggplot for single feature", {
  skip_if_not_installed("Seurat")
  skip("Requires spatial Seurat object with tissue coordinates")
  spat <- GetVisiumData()

  p <- ScSpatialFeaturePlot(spat, features = "MLANA")
  expect_s3_class(p, "ggplot")
})

test_that("ScSpatialFeaturePlot returns ggplot for multiple features", {
  skip_if_not_installed("Seurat")
  skip("Requires spatial Seurat object with tissue coordinates")
  spat <- GetVisiumData()

  p <- ScSpatialFeaturePlot(spat, features = c("MLANA", "MITF"))
  expect_s3_class(p, "ggplot")
})

test_that("combine=FALSE returns named list of correct length", {
  skip_if_not_installed("Seurat")
  skip("Requires spatial Seurat object with tissue coordinates")
  spat <- GetVisiumData()

  features <- c("MLANA", "MITF", "SOX10")
  plots <- ScSpatialFeaturePlot(spat, features = features, combine = FALSE)

  expect_type(plots, "list")
  expect_length(plots, 3)
  expect_named(plots, features)
  for (p in plots) {
    expect_s3_class(p, "ggplot")
  }
})

test_that("log transformation works", {
  skip_if_not_installed("Seurat")
  skip("Requires spatial Seurat object with tissue coordinates")
  spat <- GetVisiumData()

  p <- ScSpatialFeaturePlot(spat, features = "MLANA", logs = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("flip argument produces valid plots for all modes", {
  skip_if_not_installed("Seurat")
  skip("Requires spatial Seurat object with tissue coordinates")
  spat <- GetVisiumData()

  for (flip_val in 1:3) {
    p <- ScSpatialFeaturePlot(spat, features = "MLANA", flip = flip_val)
    expect_s3_class(p, "ggplot")
  }
})

test_that("input Seurat object is not mutated", {
  skip_if_not_installed("Seurat")
  skip("Requires spatial Seurat object with tissue coordinates")
  spat <- GetVisiumData()

  meta_before <- colnames(spat[[]])
  ScSpatialFeaturePlot(spat, features = "MLANA")
  meta_after <- colnames(spat[[]])

  expect_identical(meta_before, meta_after)
})
