test_that("GetVisiumData returns a Seurat object when data is available", {
  skip_if_not_installed("Seurat")

  rda_path <- system.file("extdata", "VisiumData.rda",
                          package = "ScSpatialFeaturePlot")
  skip_if(rda_path == "", "VisiumData.rda not bundled in test environment")

  spat <- GetVisiumData()
  expect_true(inherits(spat, "Seurat"))
})

test_that("GetVisiumData gives informative error when data is missing", {
  skip_if_not_installed("Seurat")

  rda_path <- system.file("extdata", "VisiumData.rda",
                          package = "ScSpatialFeaturePlot")
  # Only test this if the data truly isn't installed
  skip_if(rda_path != "", "VisiumData.rda is present; cannot test missing-file path")

  expect_error(
    GetVisiumData(),
    "VisiumData.rda not found"
  )
})
