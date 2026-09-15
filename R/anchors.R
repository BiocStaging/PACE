#' Contamination anchor genes per cell type
#'
#' Returns, for each cell type, the genes used as negative-control anchors when
#' identifying the per-cell contamination loading: genes owned by another cell
#' type (high mean expression there) that are near-absent in the focal type.
#' Anchors are recomputed with the same selection as the fit, and ranked by
#' expression in their owning cell type. This is the source of the manuscript's
#' supplementary anchor-gene tables.
#'
#' @param object A [PACEFit] (supplies the cell-type and image columns used at
#'   fit time).
#' @param spe The [SpatialExperiment::SpatialExperiment] that was fitted.
#' @param top_n Number of top anchors to list per cell type (default 10).
#' @param assay_name Counts assay name (default `"counts"`).
#' @return A data frame with one row per cell type: `cellType`, `n_anchors`
#'   (total anchors), and `top_anchors` (the `top_n` strongest, comma-separated).
#' @examples
#' spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' anchorGenes(fit, spe)
#' @export
anchorGenes <- function(object, spe, top_n = 10L, assay_name = "counts") {
  stopifnot(methods::is(object, "PACEFit"))
  ct     <- as.character(SummarizedExperiment::colData(spe)[[object@params$celltype_col]])
  types  <- sort(unique(ct))
  coords <- as.matrix(SpatialExperiment::spatialCoords(spe))
  Y      <- .pace_as_dgc(Matrix::t(SummarizedExperiment::assay(spe, assay_name)))   # cells x genes, sparse
  genes  <- colnames(Y)

  ## image grouping: within-image homotypic cores (constant for a single section)
  img_col <- object@params$image_col
  image <- if (!is.null(img_col) && img_col %in% colnames(SummarizedExperiment::colData(spe)))
             as.character(SummarizedExperiment::colData(spe)[[img_col]])
           else rep("all", nrow(Y))

  anchors <- pace_anchors(coords, Y, ct, image, types)
  mask    <- anchors$mask                                    # types x genes (0/1)

  ## owner mean = max mean expression of each gene across cell types
  ## (per-type colMeans computed in C++ from the sparse counts)
  type_means <- t(pace_group_column_means_cpp(Y, .pace_codes(ct, types), length(types),
                                              detection = FALSE, n_threads = 1L))  # genes x types
  dimnames(type_means) <- list(genes, types)
  owner_mean <- apply(type_means, 1, max)
  names(owner_mean) <- genes

  rows <- lapply(types, function(focal) {
    ag     <- genes[mask[focal, ] == 1]
    ranked <- ag[order(owner_mean[ag], decreasing = TRUE)]
    data.frame(cellType    = focal,
               n_anchors   = length(ag),
               top_anchors = paste(utils::head(ranked, top_n), collapse = ", "),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
