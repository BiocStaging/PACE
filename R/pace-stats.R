## pace-stats.R -- per-cell-type statistics for the variance decomposition and the
## driver scores, so paceDecompose() and paceDrivers() never hold an n x G matrix.
##
## What the reporting layers need from the fitted means, for cell type c (its n_c
## cells) and gene g:
##   mu_mean[c, g]  = mean_{i in c} mu_ig                          (V_disp, V_resid)
##   toff_var[c, g] = var_{i in c} log1p(mu_spill_ig / mu_bio_ig)   (V_spill)
##   toff_any_nonzero = any technical offset != 0                   (use_bleed gate)
## All other decomposition inputs come from the counts and the stored design:
##   per-type count means (focal_mean, spec), and per-type covariances of the
##   kernel columns, turning V_state into quadratic forms (engine-decomp.R).
##
## Sources, in order: the statistics the solver stored at fit time (fit$stats);
## the n x G matrices of an older fit that retains them; or, for a stripped older
## fit, a rebuild of mu a block of genes at a time from the fit and the counts.

## The statistics record shared by both solvers and the legacy paths.
.pace_statistics_record <- function(cell_types, n_type, mu_mean, toff_var, toff_any_nonzero) {
  list(version = 1L,
       cell_types = as.character(cell_types),
       n_type = stats::setNames(as.integer(n_type), cell_types),
       mu_mean = mu_mean,
       toff_var = toff_var,
       toff_any_nonzero = isTRUE(toff_any_nonzero))
}

## Statistics from full n x G matrices (a solver's own, or an older fit's stored
## ones). Gene by gene, so no further n x G copy is made:
##   mu_mean  = mean(mu[cells, g], na.rm = TRUE)
##   toff_var = stats::var(toff[cells, g], na.rm = TRUE)
## `celltype` is the per-cell label, `cell_types` the rows of the result.
.pace_statistics_from_matrices <- function(mu, toff, celltype, cell_types, gene_names = colnames(mu)) {
  n_genes <- ncol(mu)
  cells_by_type <- lapply(cell_types, function(type) which(celltype == type))
  mu_mean <- matrix(0, length(cell_types), n_genes, dimnames = list(cell_types, gene_names))
  toff_var <- matrix(NA_real_, length(cell_types), n_genes, dimnames = list(cell_types, gene_names))
  toff_any_nonzero <- FALSE
  for (g in seq_len(n_genes)) {
    mu_column <- mu[, g]
    toff_column <- if (is.null(toff)) numeric(nrow(mu)) else toff[, g]
    toff_any_nonzero <- toff_any_nonzero || any(toff_column != 0, na.rm = TRUE)
    for (ci in seq_along(cell_types)) {
      rows <- cells_by_type[[ci]]
      if (!length(rows)) next
      mu_mean[ci, g] <- mean(mu_column[rows], na.rm = TRUE)
      toff_var[ci, g] <- stats::var(toff_column[rows], na.rm = TRUE)
    }
  }
  .pace_statistics_record(cell_types, lengths(cells_by_type), mu_mean, toff_var, toff_any_nonzero)
}

## Statistics for a stripped fit, rebuilding mu from the fit and the counts a
## block of genes at a time (all cells, `genes_per_block` genes), so memory is
## n x genes_per_block rather than n x G. Same rebuild formula as the solver's
## final pass (.pace_mu_block()).
.pace_statistics_by_rebuild <- function(object, spe, max_block_entries = 5e6) {
  df <- object@context$df
  genes <- object@context$genes
  cell_types <- object@cellTypes
  celltype <- as.character(df$celltype)
  n_cells <- nrow(df)
  inputs <- .pace_mu_inputs(object, spe)
  genes_per_block <- max(1L, as.integer(floor(max_block_entries / max(n_cells, 1L))))
  cells_by_type <- lapply(cell_types, function(type) which(celltype == type))
  mu_mean <- matrix(0, length(cell_types), length(genes), dimnames = list(cell_types, genes))
  toff_var <- matrix(NA_real_, length(cell_types), length(genes), dimnames = list(cell_types, genes))
  toff_any_nonzero <- FALSE
  for (start in seq(1L, length(genes), by = genes_per_block)) {
    gene_idx <- start:min(start + genes_per_block - 1L, length(genes))
    block <- .pace_mu_block(object, NULL, inputs, gene_idx = gene_idx)
    if (is.null(block$mu_spill)) {
      mu_block <- block$mu_bio
      toff_block <- matrix(0, nrow(mu_block), ncol(mu_block))
    } else {
      mu_block <- pmax(block$mu_bio + block$mu_spill, 1e-6)
      toff_block <- log1p(block$mu_spill / block$mu_bio)
    }
    toff_any_nonzero <- toff_any_nonzero || any(toff_block != 0, na.rm = TRUE)
    for (ci in seq_along(cell_types)) {
      rows <- cells_by_type[[ci]]
      if (!length(rows)) next
      mu_mean[ci, gene_idx] <- apply(mu_block[rows, , drop = FALSE], 2, mean, na.rm = TRUE)
      toff_var[ci, gene_idx] <- apply(toff_block[rows, , drop = FALSE], 2, stats::var, na.rm = TRUE)
    }
    rm(block, mu_block, toff_block)
  }
  .pace_statistics_record(cell_types, lengths(cells_by_type), mu_mean, toff_var, toff_any_nonzero)
}

## The per-cell-type statistics of a fit, from the best available source.
.pace_fit_statistics <- function(object, spe = NULL) {
  f <- object@fit
  if (!is.null(f$stats)) return(f$stats)
  if (!is.null(f$mu)) {
    toff <- f$technical_offset_mat
    return(.pace_statistics_from_matrices(f$mu, toff, as.character(object@context$df$celltype),
                                          object@cellTypes, gene_names = object@context$genes))
  }
  if (is.null(spe))
    stop("this fit was saved without its fitted means or statistics; pass the ",
         "SpatialExperiment it was fitted on.", call. = FALSE)
  .pace_statistics_by_rebuild(object, spe)
}

## The fitted genes' counts as a sparse cells x genes matrix, checked against the fit.
.pace_counts_for_fit <- function(object, spe) {
  counts <- SummarizedExperiment::assay(spe, object@params$assay_name)
  df <- object@context$df
  if (ncol(counts) != nrow(df))
    stop("`spe` has ", ncol(counts), " cells but the fit has ", nrow(df),
         "; pass the same object used for paceModel().", call. = FALSE)
  .pace_as_dgc(Matrix::t(counts[object@context$genes, , drop = FALSE]))
}

## Per-group covariance matrices of the columns of `values` (n x p), computed in
## C++ exactly as stats::cov(values[group == g, , drop = FALSE]). Returns a list
## of p x p matrices named by `groups`.
.pace_group_covariances <- function(values, group, groups, threads = 1L) {
  values <- as.matrix(values)
  storage.mode(values) <- "double"
  codes <- .pace_codes(as.character(group), groups)
  covariance <- pace_group_covariances_cpp(values, codes, length(groups),
                                           n_threads = .pace_thread_count(threads))
  out <- lapply(seq_along(groups), function(g) {
    matrix(covariance[, , g], ncol(values), ncol(values),
           dimnames = list(colnames(values), colnames(values)))
  })
  stats::setNames(out, groups)
}
