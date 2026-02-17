#' Load Example Visium Spatial Transcriptomics Data
#'
#' Loads a pre-processed Seurat object containing 10X Visium spatial
#' transcriptomics data, for use in tutorials and examples.
#'
#' @return A Seurat object containing spatial transcriptomics data.
#'
#' @examples
#' \dontrun{
#' spat <- GetVisiumData()
#' }
#'
#' @export
GetVisiumData <- function() {
  rda_path <- system.file("extdata", "VisiumData.rda",
                          package = "ScSpatialFeaturePlot")
  if (rda_path == "") {
    stop("VisiumData.rda not found. The example data may not be installed.",
         call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  load(rda_path, envir = env)
  return(env$spat)
}
