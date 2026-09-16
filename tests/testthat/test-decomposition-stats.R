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
  reference <- reference_pace_top_drivers(fit@fit, shrunk@neighbourSlopes, pair$reference$blocks,
                                          fit@cellTypes, mu_means, resp_term = fit@params$resp_term,
                                          resp_dummy = fit@context$df$.resp_dummy)
  expect_tables_close(drivers, reference, label = "condition drivers")

  # the Pratt pair attribution behind pairVariance() and plotPairHeatmap()
  mv <- PACE:::.pace_mv_shim(shrunk)
  for (prefix in list(NULL, fit@params$resp_term)) {
    expect_tables_close(PACE:::pace_pair_variance_pratt(mv, cond_prefix = prefix),
                        reference_pair_variance_pratt(mv, cond_prefix = prefix),
                        label = paste("pratt", if (is.null(prefix)) "spatial" else prefix))
  }
})

test_that("group covariances match stats::cov and do not depend on threads", {
  set.seed(3)
  values <- matrix(stats::rnorm(6000), 1000, 6)
  values[, 6] <- 0                     # a constant column, as after the sparse-pair drop
  group <- sample(c("a", "b", "c"), 1000, replace = TRUE)
  groups <- c("a", "b", "c", "empty")
  one <- PACE:::.pace_group_covariances(values, group, groups, threads = 1L)
  four <- PACE:::.pace_group_covariances(values, group, groups, threads = 4L)
  expect_identical(one, four)
  expect_identical(dim(one), c(6L, 6L, 4L))
  for (g in seq_len(3)) {
    expect_tables_close(one[, , g], stats::cov(values[group == groups[g], , drop = FALSE]),
                        label = groups[g])
  }
  expect_true(all(is.nan(one[, , 4])))
})

test_that("per-group moments of a dense block match mean() and stats::var()", {
  set.seed(11)
  values <- matrix(stats::rnorm(800), 100, 8)
  values[1, 1] <- NA_real_                       # na.rm = TRUE drops it
  group <- c(rep("a", 60), rep("b", 39), "c")    # c has one cell: var is NA
  codes <- PACE:::.pace_codes(group, c("a", "b", "c", "empty"))
  one <- PACE:::pace_dense_group_moments_cpp(values, 100L, 8L, codes, 4L, TRUE, TRUE, 1L)
  four <- PACE:::pace_dense_group_moments_cpp(values, 100L, 8L, codes, 4L, TRUE, TRUE, 4L)
  expect_identical(one, four)
  expect_true(one$any_nonzero)
  for (g in seq_len(3)) {
    rows <- codes == g - 1L
    expect_tables_close(one$mean[g, ], apply(values[rows, , drop = FALSE], 2, mean, na.rm = TRUE),
                        label = "mean")
    expect_tables_close(one$variance[g, ],
                        apply(values[rows, , drop = FALSE], 2, stats::var, na.rm = TRUE),
                        label = "variance")
  }
  expect_true(all(is.na(one$variance[3, ])))     # one cell
  expect_true(all(is.na(one$variance[4, ])))     # no cells
  zero <- PACE:::pace_dense_group_moments_cpp(matrix(0, 0, 0), 100L, 8L, codes, 4L, TRUE, TRUE, 1L)
  expect_false(zero$any_nonzero)
  expect_true(all(zero$variance[1:2, ] == 0))
})

test_that("the solver's final-pass accumulation matches the R expressions", {
  set.seed(5)
  n <- 40L
  genes <- 6L
  eta <- matrix(stats::rnorm(n * genes, -2), n, genes)
  ambient <- matrix(abs(stats::rnorm(n * genes)), n, genes)
  offset <- log(stats::runif(n, 500, 2000))
  rho <- c(0, stats::runif(n - 1L, 0, 0.3))
  codes <- PACE:::.pace_codes(rep(c("a", "b", "c"), length.out = n), c("a", "b", "c"))
  pass <- PACE:::pace_final_pass_statistics_cpp(eta, offset, ambient, rho, codes, 3L, TRUE)
  mu_bio <- pmax(exp(eta + offset), 1e-6)
  mu_spill <- pmax(ambient * rho, 0)
  mu <- pmax(mu_bio + mu_spill, 1e-6)
  toff <- log1p(mu_spill / mu_bio)
  expect_identical(pass$mu, mu)
  expect_identical(pass$toff, toff)
  expect_true(pass$any_nonzero)
  expect_tables_close(pass$mu_column_sum, colSums(mu), label = "column sums")
  expect_tables_close(pass$spill_row_sum, rowSums(mu_spill), label = "spill row sums")
  expect_tables_close(pass$total_row_sum, rowSums(mu), label = "total row sums")
  for (g in seq_len(3)) {
    rows <- codes == g - 1L
    expect_tables_close(pass$mu_group_sum[g, ], colSums(mu[rows, , drop = FALSE]), label = "group sums")
    expect_tables_close(pass$toff_variance[g, ],
                        apply(toff[rows, , drop = FALSE], 2, stats::var), label = "offset variance")
  }
  expect_error(PACE:::pace_final_pass_statistics_cpp(replace(eta, 1, NA_real_), offset, ambient, rho,
                                                     codes, 3L, FALSE), "finite")
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
