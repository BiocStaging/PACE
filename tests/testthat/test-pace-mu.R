# The fitted means and the contamination log-offset are rebuilt from the fit
# rather than stored. The rebuild must be exact: it feeds the decomposition, and
# a near-miss would move published variance shares without failing.

spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))

test_that("the packaged fit re-decomposes to exactly its stored table", {
  # The packaged fit carries no n x G matrices, so paceDecompose() must rebuild
  # them. Its stored decomposition was computed at fit time from the solver's
  # own matrices, which makes this an end-to-end equality check.
  skip_if_not(is.null(fit@fit$mu), "packaged fit unexpectedly retains mu")
  stored <- fit@varianceDecomposition$perGene
  redone <- paceDecompose(fit, spe)@varianceDecomposition$perGene

  num <- vapply(stored, is.numeric, logical(1))
  expect_identical(dim(redone), dim(stored))
  # not bit-for-bit: the rebuild's BLAS products differ in the last bits across
  # platforms
  expect_equal(as.matrix(redone[num]), as.matrix(stored[num]), tolerance = 1e-10)
})

test_that("the spillover block survives the rebuild", {
  # technical_offset_mat gates the spillover block through `use_bleed`. Rebuild
  # mu but not the offset and the block silently reports zero rather than
  # erroring, so assert it is actually populated.
  sp <- paceDecompose(fit, spe)@varianceDecomposition$perGene[["Spillover %"]]
  expect_true(any(sp > 0, na.rm = TRUE))
  expect_gt(median(sp, na.rm = TRUE), 0)
})

test_that("stored matrices are returned untouched when present", {
  object <- fit
  n <- nrow(object@context$df)
  g <- length(object@context$genes)
  object@fit$mu                   <- matrix(1.5, n, g)
  object@fit$technical_offset_mat <- matrix(0.25, n, g)

  parts <- PACE:::.pace_mu_parts(object, spe)
  expect_identical(parts$mu, object@fit$mu)
  expect_identical(parts$technical_offset_mat, object@fit$technical_offset_mat)
})

test_that("the rebuild uses the fit's own settings, not the exported defaults", {
  # ambientField() defaults to h_tech = 5; a fit made with another bandwidth
  # must not be handed that field. Perturbing the recorded value must move the
  # result, or the fit's settings are being ignored.
  base <- PACE:::.pace_mu_parts(fit, spe)
  other <- fit
  other@params$h_tech <- 12
  moved <- PACE:::.pace_mu_parts(other, spe)
  expect_false(isTRUE(all.equal(base$mu, moved$mu)))

  other2 <- fit
  other2@params$edge_correct <- !fit@params$edge_correct
  expect_false(isTRUE(all.equal(base$mu, PACE:::.pace_mu_parts(other2, spe)$mu)))
})

test_that("a fit predating edge_correct is refused rather than guessed at", {
  object <- fit
  object@params$edge_correct <- NULL
  expect_error(PACE:::.pace_mu_parts(object, spe), "edge_correct")
})

test_that("contamination = none rebuilds no spillover", {
  object <- fit
  object@params$contamination   <- "none"
  object@fit$percell_bleed_rho  <- NULL
  parts <- PACE:::.pace_mu_parts(object, spe)
  expect_true(all(parts$technical_offset_mat == 0))
  expect_true(all(parts$mu >= 1e-6))
})

test_that("a mismatched spe is refused", {
  expect_error(PACE:::.pace_mu_parts(fit, spe[, 1:100]), "cells but the fit has")
})

test_that("the packaged fit re-scores to exactly its stored driver tables", {
  # The driver scores read the fitted means too. Without a rebuild a stripped
  # fit returned an empty list behind a warning.
  skip_if_not(is.null(fit@fit$mu), "packaged fit unexpectedly retains mu")
  # a failure inside the scoring still surfaces as a warning and an empty list
  expect_no_warning(redone <- paceDrivers(fit, spe)@topDrivers)
  expect_gt(length(fit@topDrivers), 0)
  expect_gt(sum(vapply(redone, function(d) nrow(d$scores), integer(1))), 0)
  expect_identical(names(redone), names(fit@topDrivers))
  for (k in names(fit@topDrivers)) {
    # not bit-for-bit: the rebuild's BLAS products differ in the last bits
    # across platforms
    expect_equal(redone[[k]]$scores, fit@topDrivers[[k]]$scores, tolerance = 1e-10)
    expect_identical(redone[[k]]$status, fit@topDrivers[[k]]$status)
  }
})

test_that("a stripped fit without spe is refused, not scored empty", {
  skip_if_not(is.null(fit@fit$mu), "packaged fit unexpectedly retains mu")
  expect_error(paceDrivers(fit), "pass the SpatialExperiment")
  # pairs given positionally used to land in `spe`
  expect_error(paceDrivers(fit, list(c("Macrophage", "Tumour"))), "pass pairs by name")
})

test_that("chunked cell-type means match the means of the whole rebuilt mu", {
  # paceDrivers() averages mu a block of cells at a time so the n x G matrix is
  # never held. Block boundaries must not change the answer, including a block
  # size that does not divide the number of cells.
  full_mu <- PACE:::.pace_mu_parts(fit, spe)$mu
  types <- fit@cellTypes
  intercepts <- fit@fit$re_meta$Z[, paste0(types, "::(Intercept)"), drop = FALSE]
  reference <- t(vapply(types, function(ct)
    colMeans(full_mu[intercepts[, paste0(ct, "::(Intercept)")] != 0, , drop = FALSE]),
    numeric(ncol(full_mu))))

  for (chunk in c(nrow(full_mu), 1000L, 777L)) {
    means <- PACE:::.pace_mu_means_by_celltype(fit, spe, chunk_size = chunk)
    expect_identical(dimnames(means), dimnames(reference))
    expect_equal(means, reference, tolerance = 1e-10)
  }

  # a fit that stores mu is averaged directly, to the same values
  stored <- fit
  stored@fit$mu <- full_mu
  expect_equal(PACE:::.pace_mu_means_by_celltype(stored), reference, tolerance = 1e-10)
})
