# The variance decomposition, single-frame table and driver scores computed from
# per-cell-type statistics must reproduce the per-cell implementation
# (helper-decomposition-reference.R) on single-section, condition and
# contamination = "none" fits, and the three statistics sources (stored at fit
# time, stored matrices, stripped rebuild) must agree.

small_bc <- function(fraction = 0.25) {
  f <- system.file("extdata", "bc_xenium_subset.rds", package = "PACE")
  skip_if(f == "", "example dataset not installed")
  spe <- readRDS(f)
  x <- SpatialExperiment::spatialCoords(spe)[, 1]
  spe[, x <= stats::quantile(x, fraction)]
}

small_mel <- function() {
  f <- system.file("extdata", "mel_cosmx_subset.rds", package = "PACE")
  skip_if(f == "", "example dataset not installed")
  spe <- readRDS(f)
  spe[, spe$image %in% c("32134_26", "32151_13", "32158_19")]
}

# Decomposition tables through paceDecompose() for a fit, and through the
# per-cell reference on the same fit's stored matrices.
decomposition_pair <- function(fit, spe) {
  new <- paceDecompose(fit, spe)@varianceDecomposition
  df <- fit@context$df
  Y <- t(as.matrix(SummarizedExperiment::assay(spe, "counts")))[, fit@context$genes, drop = FALSE]
  blocks <- reference_pace_decompose(fit@fit, df, Y, fit@cellTypes, fit@context$X_fixed,
                                     resp_term = fit@params$resp_term)
  block <- if (is.null(fit@params$condition_col)) blocks$gene_focal_4block else blocks$gene_focal_5block
  per_gene <- reference_single_frame_decomp_obs(Y, df[[fit@params$celltype_col]],
                                                as.numeric(Matrix::rowSums(Y)), block)
  list(new = new, reference = list(perGene = per_gene, blocks = blocks))
}

check_fit <- function(fit, spe, label) {
  expect_false(is.null(fit@fit$stats), label = label)
  pair <- decomposition_pair(fit, spe)
  expect_tables_close(pair$new, pair$reference, label = paste(label, "stats vs per-cell"))

  legacy_matrices <- fit
  legacy_matrices@fit$stats <- NULL
  expect_tables_close(paceDecompose(legacy_matrices, spe)@varianceDecomposition, pair$new,
                      label = paste(label, "stored matrices"))

  stripped <- fit
  stripped@fit[c("stats", "mu", "technical_offset_mat", "bleed_offset_mat")] <- NULL
  expect_tables_close(paceDecompose(stripped, spe)@varianceDecomposition, pair$new,
                      label = paste(label, "stripped rebuild"))
  invisible(pair)
}

test_that("single-section fits decompose identically from statistics", {
  skip_if_not_installed("SpatialExperiment")
  spe <- small_bc()
  # the data-informed prior warns on the tiny crop's absent kernel columns
  fit <- suppressWarnings(paceModel(spe, celltype_col = "cellType", n_iter = 3L, threads = 2L,
                                    verbose = FALSE))
  check_fit(fit, spe, "bc")
})

test_that("contamination = 'none' fits decompose identically from statistics", {
  skip_if_not_installed("SpatialExperiment")
  spe <- small_bc()
  fit <- suppressWarnings(paceModel(spe, celltype_col = "cellType", contamination = "none",
                                    n_iter = 2L, threads = 2L, verbose = FALSE))
  pair <- check_fit(fit, spe, "none")
  expect_false(fit@fit$stats$toff_any_nonzero)
  expect_true(any(pair$new$blocks$gene_focal_5block$V_spill > 0))
})

test_that("condition fits with an image block decompose and score identically", {
  skip_if_not_installed("SpatialExperiment")
  spe <- small_mel()
  fit <- paceModel(spe, celltype_col = "cellType", condition_col = "Responder",
                   image_col = "image", kernel_per_image = TRUE, image_re = "intercept",
                   n_iter = 3L, threads = 2L, verbose = FALSE)
  pair <- check_fit(fit, spe, "condition")
  expect_true(any(pair$new$blocks$gene_focal_5block$V_state_responder > 0))

  shrunk <- paceShrink(fit, null_correlation = FALSE)
  shrunk@varianceDecomposition <- pair$new
  drivers <- paceDrivers(shrunk, spe)@topDrivers
  # reference: the driver scores with per-type means of the stored mu
  cells_by_type <- lapply(fit@cellTypes, function(type) which(as.character(fit@context$df$celltype) == type))
  mu_means <- t(vapply(cells_by_type, function(rows) colMeans(fit@fit$mu[rows, , drop = FALSE]),
                       numeric(ncol(fit@fit$mu))))
  dimnames(mu_means) <- list(fit@cellTypes, fit@context$genes)
  reference <- PACE:::pace_top_drivers(fit@fit, shrunk@neighbourSlopes, pair$reference$blocks,
                                       fit@cellTypes, mu_means, resp_term = fit@params$resp_term,
                                       resp_dummy = fit@context$df$.resp_dummy)
  expect_tables_close(drivers, reference, label = "condition drivers")
})

test_that("group covariances match stats::cov and do not depend on threads", {
  set.seed(3)
  values <- matrix(stats::rnorm(6000), 1000, 6)
  values[, 6] <- 0                     # a constant column, as after the sparse-pair drop
  group <- sample(c("a", "b", "c"), 1000, replace = TRUE)
  one <- PACE:::.pace_group_covariances(values, group, c("a", "b", "c", "empty"), threads = 1L)
  four <- PACE:::.pace_group_covariances(values, group, c("a", "b", "c", "empty"), threads = 4L)
  expect_identical(one, four)
  for (g in c("a", "b", "c")) {
    expect_tables_close(one[[g]], stats::cov(values[group == g, , drop = FALSE]), label = g)
  }
  expect_true(all(is.nan(one$empty)))
})

test_that("single-frame statistics match the dense computation and do not depend on threads", {
  skip_if_not_installed("SpatialExperiment")
  spe <- small_bc()
  Y <- t(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  counts <- PACE:::.pace_as_dgc(Y)
  n_count <- as.numeric(Matrix::rowSums(counts))
  types <- sort(unique(spe$cellType))
  codes <- PACE:::.pace_codes(as.character(spe$cellType), types)
  one <- PACE:::pace_single_frame_statistics_cpp(counts, n_count, codes, length(types), 1L)
  four <- PACE:::pace_single_frame_statistics_cpp(counts, n_count, codes, length(types), 4L)
  expect_identical(one, four)
  y_log <- log1p(Y * (1e4 / n_count))
  for (ti in seq_along(types)) {
    focal <- y_log[codes == ti - 1L, , drop = FALSE]
    focal_mean <- colMeans(focal)
    expect_tables_close(one$focal_mean[ti, ], focal_mean, label = types[ti])
    expect_tables_close(one$within_ss[ti, ], colSums(sweep(focal, 2, focal_mean, "-")^2), label = types[ti])
  }
  expect_tables_close(one$global_mean, colMeans(y_log), label = "global")
  expect_error(PACE:::pace_single_frame_statistics_cpp(counts, replace(n_count, 1, 0), codes,
                                                       length(types), 1L), "library size")
})

test_that("cells with zero total counts are refused at fit time", {
  skip_if_not_installed("SpatialExperiment")
  spe <- small_bc(0.1)
  counts <- SummarizedExperiment::assay(spe, "counts")
  counts[, 1] <- 0
  SummarizedExperiment::assay(spe, "counts") <- counts
  expect_error(paceModel(spe, celltype_col = "cellType", n_iter = 1L, threads = 1L, verbose = FALSE),
               "zero total counts")
})
