#' Spatial Feature Plot with Shared Color Scale
#'
#' An extension of Seurat's \code{SpatialFeaturePlot} that allows plotting
#' multiple genes on the same color scale, facilitating comparison of expression
#' magnitude across genes. Also supports rotation of tissue coordinates.
#'
#' @param spat A Seurat object containing spatial transcriptomics data.
#' @param features Character vector of gene names or metadata column names to plot.
#' @param slot Character string specifying which data slot to pull from.
#'   Default is \code{"counts"}.
#' @param facet_labeller A labeller function for facet labels. Default is \code{NULL}.
#' @param logs Logical; if \code{TRUE}, log-transform expression values using
#'   \code{log(x + 1)}. Default is \code{FALSE}.
#' @param pt.size.factor Numeric; point size multiplier for the scatter plot.
#'   Default is \code{1.6}.
#' @param cells Character vector of cell names to include in the plot. Default
#'   is all cells in \code{spat}.
#' @param LegLabel Character string for the legend title. Default is \code{"Count"}.
#' @param uq Numeric; upper quantile for value clipping. Default is \code{0.999}.
#' @param lq Numeric; lower quantile for value clipping. Default is \code{0}.
#' @param flip Integer specifying rotation mode: 1 = flip x-axis, 2 = transpose
#'   axes, 3 = flip y-axis. Default is \code{3}.
#' @param ncol Integer; number of columns in the faceted plot. Default is
#'   \code{NULL} (automatic).
#' @param nrow Integer; number of rows in the faceted plot. Default is
#'   \code{NULL} (automatic).
#' @param combine Logical; if \code{TRUE}, return a single combined faceted plot.
#'   If \code{FALSE}, return a list of individual ggplot objects. Default is \code{TRUE}.
#'
#' @return A \code{ggplot} object (if \code{combine = TRUE}) or a list of
#'   \code{ggplot} objects (if \code{combine = FALSE}).
#'
#' @examples
#' \dontrun{
#' # Single gene
#' ScSpatialFeaturePlot(spat, features = "MLANA", pt.size.factor = 2)
#'
#' # Multiple genes on the same color scale
#' ScSpatialFeaturePlot(spat, features = c("MLANA", "MITF"), pt.size.factor = 2)
#'
#' # Independent scales per gene
#' plots <- ScSpatialFeaturePlot(spat, features = c("MLANA", "MITF"), combine = FALSE)
#' do.call(gridExtra::grid.arrange, c(plots, ncol = 2))
#' }
#'
#' @importFrom Seurat GetTissueCoordinates FetchData Cells MinMax
#' @importFrom ggplot2 ggplot geom_point aes scale_color_viridis_c theme_bw
#'   theme element_blank labs facet_wrap
#' @importFrom rlang .data
#' @export
ScSpatialFeaturePlot <- function(spat, features, slot = "counts",
                                  facet_labeller = NULL,
                                  logs = FALSE, pt.size.factor = 1.6,
                                  cells = Seurat::Cells(spat),
                                  LegLabel = "Count",
                                  uq = 0.999, lq = 0,
                                  flip = 3, ncol = NULL, nrow = NULL,
                                  combine = TRUE) {

  # --- Input validation ---

  if (!inherits(spat, "Seurat")) {
    stop("'spat' must be a Seurat object.", call. = FALSE)
  }
  if (!is.character(features) || length(features) == 0) {
    stop("'features' must be a non-empty character vector.", call. = FALSE)
  }
  if (!flip %in% c(1, 2, 3)) {
    stop("'flip' must be 1, 2, or 3.", call. = FALSE)
  }
  if (!is.numeric(uq) || !is.numeric(lq) || uq < lq) {
    stop("'uq' and 'lq' must be numeric with uq >= lq.", call. = FALSE)
  }

  # --- Handle combine = FALSE: return list of individual plots ---
  if (!combine) {
    plotlist <- lapply(features, function(one.feature) {
      ScSpatialFeaturePlot(spat, one.feature, slot, facet_labeller,
                           logs, pt.size.factor, cells, LegLabel,
                           uq, lq, flip, ncol = 1, nrow = 1,
                           combine = TRUE) +
        ggplot2::theme(legend.position = "top")
    })
    names(plotlist) <- features
    return(plotlist)
  }

  # --- Quantile clipping helper ---
  MinMaxq <- function(X, uq = 0.99, lq = 0.05) {
    Seurat::MinMax(X,
                   min = quantile(X, lq, na.rm = TRUE),
                   max = quantile(X, uq, na.rm = TRUE))
  }

  # --- Get and transform tissue coordinates ---
  rawcoords <- Seurat::GetTissueCoordinates(spat)[cells, ]
  if (flip == 1) {
    rawcoords$x <- -rawcoords$imagecol
    rawcoords$y <- rawcoords$imagerow
  } else if (flip == 2) {
    rawcoords$x <- rawcoords$imagerow
    rawcoords$y <- rawcoords$imagecol
  } else {
    rawcoords$x <- rawcoords$imagecol
    rawcoords$y <- -rawcoords$imagerow
  }

  # --- Build plotting data frame (no mutation of input object) ---
  meta_cols <- colnames(spat[[]])
  plot_dfs <- lapply(features, function(fname) {
    if (fname %in% meta_cols) {
      feat_values <- spat[[fname]][cells, 1]
    } else {
      feat_values <- Seurat::FetchData(spat, layer = slot, vars = fname)[cells, 1]
    }
    feat_clipped <- MinMaxq(feat_values, uq = uq, lq = lq)

    data.frame(
      x = rawcoords$x,
      y = rawcoords$y,
      feature = feat_clipped[rownames(rawcoords)],
      fname = rep(fname, nrow(rawcoords)),
      stringsAsFactors = FALSE
    )
  })
  all.plotdf <- do.call(rbind, plot_dfs)

  # --- Optional log transform ---
  if (logs) {
    all.plotdf$feature <- log(all.plotdf$feature + 1)
  }

  # --- Build ggplot ---
  p1 <- ggplot2::ggplot(all.plotdf) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = .data$x, y = .data$y, color = .data$feature),
      size = pt.size.factor
    ) +
    ggplot2::scale_color_viridis_c(option = "turbo") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "", y = "", color = LegLabel) +
    ggplot2::facet_wrap(~fname, ncol = ncol, nrow = nrow)

  return(p1)
}
