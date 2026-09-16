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
## The two moments themselves are computed in C++ (pace::dense_group_moments);
## the functions here select the source and label the result.

## The statistics record shared by both solvers and the legacy paths.
.pace_statistics_record <- function(cell_types, n_type, mu_mean, toff_var, toff_any_nonzero) {
  list(version = 1L,
       cell_types = as.character(cell_types),
       n_type = stats::setNames(as.integer(n_type), cell_types),
       mu_mean = mu_mean,
       toff_var = toff_var,
       toff_any_nonzero = isTRUE(toff_any_nonzero))
}

## Per-cell-type moments of one dense block: the mean of `mu` and the variance of
## `toff` over each type's cells, per gene, as mean()/stats::var() computed them.
## `toff = NULL` stands for an all-zero offset and costs no allocation.
.pace_block_moments <- function(mu, toff, type_code, n_types, threads = 1L) {
  threads <- .pace_thread_count(threads)
  empty <- matrix(0, 0, 0)
  mu_stats <- pace_dense_group_moments_cpp(mu, nrow(mu), ncol(mu), type_code, n_types,
                                           want_mean = TRUE, want_variance = FALSE,
                                           n_threads = threads)
  toff_stats <- pace_dense_group_moments_cpp(if (is.null(toff)) empty else toff,
                                             nrow(mu), ncol(mu), type_code, n_types,
                                             want_mean = FALSE, want_variance = TRUE,
                                             n_threads = threads)
  list(mu_mean = mu_stats$mean, toff_var = toff_stats$variance,
       toff_any_nonzero = isTRUE(toff_stats$any_nonzero))
}

## Statistics from full n x G matrices (a solver's own, or an older fit's stored
## ones). `celltype` is the per-cell label, `cell_types` the rows of the result.
## A cell type with no cells keeps the mean 0 the fit record has always held.
.pace_statistics_from_matrices <- function(mu, toff, celltype, cell_types, gene_names = colnames(mu),
                                           threads = 1L) {
  type_code <- .pace_codes(celltype, cell_types)
  moments <- .pace_block_moments(mu, toff, type_code, length(cell_types), threads)
  n_type <- tabulate(type_code + 1L, nbins = length(cell_types))
  moments$mu_mean[n_type == 0L, ] <- 0
  dimnames(moments$mu_mean) <- list(cell_types, gene_names)
  dimnames(moments$toff_var) <- list(cell_types, gene_names)
  .pace_statistics_record(cell_types, n_type, moments$mu_mean, moments$toff_var,
                          moments$toff_any_nonzero)
}

## Statistics for a stripped fit, rebuilding mu from the fit and the counts a
## block of genes at a time (all cells, `genes_per_block` genes), so memory is
## n x genes_per_block rather than n x G. Same rebuild formula as the solver's
## final pass (.pace_mu_block()).
.pace_statistics_by_rebuild <- function(object, spe, max_block_entries = 5e6, threads = 1L) {
  df <- object@context$df
  genes <- object@context$genes
  cell_types <- object@cellTypes
  n_cells <- nrow(df)
  inputs <- .pace_mu_inputs(object, spe)
  genes_per_block <- max(1L, as.integer(floor(max_block_entries / max(n_cells, 1L))))
  type_code <- .pace_codes(as.character(df$celltype), cell_types)
  n_type <- tabulate(type_code + 1L, nbins = length(cell_types))
  mu_mean <- matrix(0, length(cell_types), length(genes), dimnames = list(cell_types, genes))
  toff_var <- matrix(NA_real_, length(cell_types), length(genes), dimnames = list(cell_types, genes))
  toff_any_nonzero <- FALSE
  for (start in seq(1L, length(genes), by = genes_per_block)) {
    gene_idx <- start:min(start + genes_per_block - 1L, length(genes))
    block <- .pace_mu_block(object, inputs, gene_idx = gene_idx)
    if (is.null(block$mu_spill)) {
      mu_block <- block$mu_bio
      toff_block <- NULL
    } else {
      mu_block <- pmax(block$mu_bio + block$mu_spill, 1e-6)
      toff_block <- log1p(block$mu_spill / block$mu_bio)
    }
    moments <- .pace_block_moments(mu_block, toff_block, type_code, length(cell_types), threads)
    mu_mean[, gene_idx] <- moments$mu_mean
    toff_var[, gene_idx] <- moments$toff_var
    toff_any_nonzero <- toff_any_nonzero || moments$toff_any_nonzero
    rm(block, mu_block, toff_block)
  }
  mu_mean[n_type == 0L, ] <- 0
  .pace_statistics_record(cell_types, n_type, mu_mean, toff_var, toff_any_nonzero)
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
  ## subset only when the fit used a subset: `counts[genes, ]` is a full copy
  if (!identical(rownames(counts), object@context$genes))
    counts <- counts[object@context$genes, , drop = FALSE]
  .pace_as_dgc(Matrix::t(counts))
}

## Per-group covariance matrices of the columns of `values` (n x p), computed in
## C++ exactly as stats::cov(values[group == g, , drop = FALSE]). Returns the
## p x p x n_groups array the decomposition core reads, slice g being group g.
.pace_group_covariances <- function(values, group, groups, threads = 1L) {
  values <- as.matrix(values)
  storage.mode(values) <- "double"
  codes <- .pace_codes(as.character(group), groups)
  pace_group_covariances_cpp(values, codes, length(groups),
                             n_threads = .pace_thread_count(threads))
}
