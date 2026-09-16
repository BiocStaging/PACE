# The decomposition and driver scores read per-cell-type statistics of the
# fitted means. A fit made by this version stores them; an older fit has them
# computed from its stored n x G matrices, or rebuilt from the fit and the counts
# when those were stripped. The packaged example fit is such a stripped older fit,
# so these tests exercise the rebuild, which must be exact: it feeds the
# decomposition, and a near-miss would move published variance shares without failing.

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

test_that("stored statistics are used as-is, and stored matrices when statistics are absent", {
  object <- fit
  sentinel <- list(version = 1L, mu_mean = "stored")
  object@fit$stats <- sentinel
  expect_identical(PACE:::.pace_fit_statistics(object, spe), sentinel)

  object <- fit
  n <- nrow(object@context$df)
  g <- length(object@context$genes)
  object@fit$mu                   <- matrix(1.5, n, g)
  object@fit$technical_offset_mat <- matrix(0.25, n, g)
  stats <- PACE:::.pace_fit_statistics(object)
  expect_true(all(stats$mu_mean == 1.5))
  expect_true(all(stats$toff_var == 0, na.rm = TRUE))
  expect_true(stats$toff_any_nonzero)
  expect_identical(rownames(stats$mu_mean), fit@cellTypes)
})

test_that("the rebuild uses the fit's own settings, not the exported defaults", {
  # ambientField() defaults to h_tech = 5; a fit made with another bandwidth
  # must not be handed that field. Perturbing the recorded value must move the
  # result, or the fit's settings are being ignored.
  base <- PACE:::.pace_fit_statistics(fit, spe)
  other <- fit
  other@params$h_tech <- 12
  moved <- PACE:::.pace_fit_statistics(other, spe)
  expect_false(isTRUE(all.equal(base$mu_mean, moved$mu_mean)))
  expect_false(isTRUE(all.equal(base$toff_var, moved$toff_var)))

  other2 <- fit
  other2@params$edge_correct <- !fit@params$edge_correct
  expect_false(isTRUE(all.equal(base$toff_var, PACE:::.pace_fit_statistics(other2, spe)$toff_var)))
})

test_that("a fit predating edge_correct is refused rather than guessed at", {
  object <- fit
  object@params$edge_correct <- NULL
  expect_error(PACE:::.pace_fit_statistics(object, spe), "edge_correct")
})

test_that("contamination = none rebuilds no spillover", {
  object <- fit
  object@params$contamination   <- "none"
  object@fit$percell_bleed_rho  <- NULL
  stats <- PACE:::.pace_fit_statistics(object, spe)
  expect_false(stats$toff_any_nonzero)
  expect_true(all(stats$toff_var == 0, na.rm = TRUE))
  expect_true(all(stats$mu_mean >= 1e-6))
})

test_that("a mismatched spe is refused", {
  expect_error(PACE:::.pace_fit_statistics(fit, spe[, 1:100]), "cells but the fit has")
  expect_error(paceDecompose(fit, spe[, 1:100]), "cells but the fit has")
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

test_that("the gene-block rebuild does not depend on the block size", {
  # A stripped fit is rebuilt a block of genes at a time, so the n x G matrix is
  # never held. Block boundaries must not change the statistics, including a
  # block that does not divide the number of genes, and they must equal the
  # statistics of the whole rebuilt matrices.
  n <- nrow(fit@context$df)
  inputs <- PACE:::.pace_mu_inputs(fit, spe)
  whole <- PACE:::.pace_mu_block(fit, inputs)
  whole_mu <- pmax(whole$mu_bio + whole$mu_spill, 1e-6)
  whole_toff <- log1p(whole$mu_spill / whole$mu_bio)
  reference <- PACE:::.pace_statistics_from_matrices(whole_mu, whole_toff,
                                                     as.character(fit@context$df$celltype),
                                                     fit@cellTypes, fit@context$genes)
  for (genes_per_block in c(length(fit@context$genes), 100L, 37L, 1L)) {
    rebuilt <- PACE:::.pace_statistics_by_rebuild(fit, spe, max_block_entries = n * genes_per_block)
    expect_identical(rebuilt$toff_any_nonzero, reference$toff_any_nonzero)
    expect_identical(dimnames(rebuilt$mu_mean), dimnames(reference$mu_mean))
    expect_equal(rebuilt$mu_mean, reference$mu_mean, tolerance = 1e-10)
    expect_equal(rebuilt$toff_var, reference$toff_var, tolerance = 1e-10)
  }
})
