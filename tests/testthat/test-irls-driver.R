# The outer PQL iteration runs in the core (pace::run_irls_loop). These tests
# pin it to a digest frozen from the build that was proven identical, field for
# field, to the R loop the core was transcribed from -- on the two canonical
# cohorts as well as on these fixtures.
#
# The oracle used to be the R loop itself, reached with R_IRLS_R_LOOP. That loop
# is gone: it was reachable only under fuse_rho, and fuse_rho existed only to
# fold the rho accumulation into the solve pass so that eta was not rebuilt a
# third time per iteration -- a cost the R loop created and the core does not
# have. Comparing against it had also become partly circular, since the R path's
# three chunk passes call the same core functions the driver does.
#
# Bit-identity is the right gate, not a tolerance: the driver makes the same
# core calls in the same order on the same buffers, so any difference at all
# means a call moved, a guard changed or a floor was dropped -- and this fit
# produced numbers in a manuscript under review.
#
# To regenerate after a DELIBERATE change, rebuild the reference with
# irls_digest() over the same five fits and say in the commit why every digest
# that moved was expected to move.

reference_fits <- function() {
  path <- testthat::test_path("irls-driver-reference.rds")
  skip_if(!file.exists(path), "the frozen reference is not installed")
  readRDS(path)
}

bc_subset <- function() {
  skip_if_not_installed("SpatialExperiment")
  path <- system.file("extdata", "bc_xenium_subset.rds", package = "PACE")
  skip_if(path == "", "example dataset not installed")
  readRDS(path)
}

test_that("the core driver reproduces the frozen fit on the breast cancer subset", {
  spe <- bc_subset()
  fit <- paceModel(spe, celltype_col = "cellType", contamination = "percell_hc",
                   dispersion = "nb1", n_iter = 6L, min_iter = 3L, threads = 2L,
                   verbose = FALSE)@fit
  expect_identical(irls_digest(fit), reference_fits()$bc)
  # A fit that stopped at iteration one would match the digest without ever
  # exercising the early stop, the dispersion freeze or the tau update.
  expect_gt(fit$n_iter, 1L)
})

test_that("the core driver reproduces the frozen fit under the other tau shrinkages", {
  spe <- bc_subset()
  reference <- reference_fits()
  # "hierarchical" is the one that transposes: R holds a block's components as a
  # K_terms x K_groups matrix filled BY ROW from the flat slice, and handing the
  # slice to the shrinkage as it lies has been a bug before.
  for (shrinkage in c("shared", "hierarchical", "half_cauchy")) {
    fit <- paceModel(spe, celltype_col = "cellType", tau_shrinkage = shrinkage,
                     n_iter = 5L, min_iter = 3L, threads = 2L, verbose = FALSE)@fit
    expect_identical(irls_digest(fit), reference[[paste0("tau_", shrinkage)]])
  }
})

test_that("the core driver reproduces the frozen fit on the Gaussian family", {
  # The identity link takes the coefficient-based convergence metric and the
  # residual-variance dispersion step, neither of which the count path reaches.
  set.seed(11)
  n <- 700L
  n_genes <- 16L
  df <- data.frame(celltype = sample(c("A", "B", "C"), n, replace = TRUE),
                   A = runif(n), B = runif(n), C = runif(n),
                   stringsAsFactors = FALSE)
  values <- matrix(rnorm(n * n_genes, mean = 4, sd = 1.5), n, n_genes)
  colnames(values) <- paste0("p", seq_len(n_genes))
  counts <- methods::as(methods::as(methods::as(values, "dMatrix"), "generalMatrix"),
                        "CsparseMatrix")
  fit <- fit_pace_mvpql_streaming(
    Y = counts, X_fixed = matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)")),
    df = df, re_specs = list(list(group_col = "celltype", formula = ~ 0 + A + B + C)),
    family = "gaussian", n_iter = 6L, min_iter = 3L, chunk_size = 8L,
    n_threads = 2L, verbose = FALSE)
  expect_identical(irls_digest(fit), reference_fits()$gaussian)
})

test_that("a binding tau cap warns from the driver, and the fit is unchanged", {
  spe <- bc_subset()
  seen <- character(0)
  fit <- NULL
  withCallingHandlers(
    fit <- paceModel(spe, celltype_col = "cellType", tau_max = 0.001, n_iter = 4L,
                     min_iter = 2L, threads = 2L, verbose = FALSE)@fit,
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  reference <- reference_fits()
  # The warnings are raised from R, after the binding has returned: raising them
  # from C++ let options(warn = 2) longjmp past the destructors that hand the
  # counts, the ambient field and Z back from R's precious list.
  expect_gt(length(seen), 0L)
  expect_identical(seen, reference$tau_cap_warnings)
  expect_identical(irls_digest(fit), reference$tau_cap)
})

test_that("the Gaussian family solves its interior in double precision", {
  # paceModelGaussian() passes interior_precision = 1, and the float per-gene
  # Cholesky is borderline on continuous intensities -- the NaN guard would fire
  # every iteration. fit_pass1 forces 0 for this family so a caller cannot lose
  # it. A float interior would show up as re-solve chatter under verbose.
  set.seed(3)
  n <- 400L
  n_genes <- 8L
  df <- data.frame(celltype = sample(c("A", "B"), n, replace = TRUE),
                   A = runif(n), B = runif(n), stringsAsFactors = FALSE)
  values <- matrix(abs(rnorm(n * n_genes, mean = 3)), n, n_genes)
  colnames(values) <- paste0("p", seq_len(n_genes))
  counts <- methods::as(methods::as(methods::as(values, "dMatrix"), "generalMatrix"),
                        "CsparseMatrix")
  chatter <- utils::capture.output(
    fit <- fit_pace_mvpql_streaming(
      Y = counts, X_fixed = matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)")),
      df = df, re_specs = list(list(group_col = "celltype", formula = ~ 0 + A + B)),
      family = "gaussian", interior_precision = 1L, n_iter = 4L, min_iter = 2L,
      n_threads = 2L, verbose = TRUE))
  expect_false(any(grepl("nan-guard", chatter, fixed = TRUE)))
  expect_true(all(is.finite(fit$U)))
})

test_that("the streamed ambient field agrees with the cached one", {
  # "stream" recomputes each chunk's ambient columns from W and Y rather than
  # slicing a cached n x G product, which is what makes the full Xenium 5K panel
  # fit in memory at 1.2M cells. The two are BIT-IDENTICAL: the chunking is
  # exact because a sparse product is column-independent, and the accumulation
  # matches CHOLMOD's once FMA contraction is disabled (fp_no_contract.hpp in
  # sparse_product.cpp -- without it clang contracts the multiply-add and the
  # answers drift by about 1e-6 in the coefficients).
  #
  # This once differed by 2.13 in U, not 1e-6, because rho_accumulate derived
  # the anchor mask's column from the AMBIENT block's gene offset. That offset
  # is the global gene index only while the block is a slice of the full panel;
  # a per-chunk block starts at zero, so every chunk after the first read the
  # wrong anchors. The offset now comes from the counts block, which is a slice
  # of the full panel in both modes.
  skip_if_not_installed("SpatialExperiment")
  path <- system.file("extdata", "bc_xenium_subset.rds", package = "PACE")
  skip_if(path == "", "example dataset not installed")
  spe <- readRDS(path)

  fit_in <- function(mode) {
    set.seed(1)
    paceModel(spe, celltype_col = "cellType", contamination = "percell_hc",
              dispersion = "nb1", n_iter = 5L, min_iter = 3L, threads = 2L,
              ambient_mode = mode, verbose = FALSE)@fit
  }
  cached <- fit_in("cache")
  streamed <- fit_in("stream")

  expect_identical(streamed$U, cached$U)
  expect_identical(streamed$B, cached$B)
  expect_identical(streamed$percell_bleed_rho, cached$percell_bleed_rho)
  expect_identical(streamed$n_iter, cached$n_iter)
})
