# The outer PQL iteration now runs in the core (pace::run_irls_loop). The R loop
# it was transcribed from is still in fit_pace_mvpql_streaming(), reachable with
# R_IRLS_R_LOOP, and it is the oracle: these tests fit the same data both ways
# and require the two to agree BIT FOR BIT, not to a tolerance.
#
# Bit-identity is the right gate here, not 1e-10. The driver makes the same core
# calls in the same order on the same buffers, so any difference at all means a
# call moved, a guard changed, or a floor was dropped -- and this fit produced
# numbers in a manuscript under review.

fit_through_r_loop <- function(fit_expression) {
  previous <- Sys.getenv("R_IRLS_R_LOOP", unset = NA_character_)
  Sys.setenv(R_IRLS_R_LOOP = "1")
  on.exit({
    if (is.na(previous)) Sys.unsetenv("R_IRLS_R_LOOP") else Sys.setenv(R_IRLS_R_LOOP = previous)
  }, add = TRUE)
  force(fit_expression)
}

# Everything the loop owns. Comparing the whole list rather than a few matrices
# is deliberate: the convergence history decides WHICH iteration is the last one,
# so a history that drifts is a fit that will drift on the next dataset.
loop_outputs <- function(fit) {
  fit[c("B", "U", "se_B", "se_U", "alpha", "tau_g_array", "tau_blocks",
        "percell_bleed_rho", "history", "n_iter", "converged")]
}

test_that("the core driver reproduces the R loop on the breast cancer subset", {
  skip_if_not_installed("SpatialExperiment")
  path <- system.file("extdata", "bc_xenium_subset.rds", package = "PACE")
  skip_if(path == "", "example dataset not installed")
  spe <- readRDS(path)

  fit_once <- function() {
    paceModel(spe, celltype_col = "cellType", contamination = "percell_hc",
              dispersion = "nb1", n_iter = 6L, min_iter = 3L, threads = 2L,
              verbose = FALSE)@fit
  }
  driven <- fit_once()
  reference <- fit_through_r_loop(fit_once())

  expect_identical(loop_outputs(driven), loop_outputs(reference))
  # A fit that stopped at iteration one would pass the comparison above without
  # exercising the early stop, the dispersion freeze or the tau update at all.
  expect_gt(driven$n_iter, 1L)
})

test_that("the core driver reproduces the R loop under the other tau shrinkages", {
  skip_if_not_installed("SpatialExperiment")
  path <- system.file("extdata", "bc_xenium_subset.rds", package = "PACE")
  skip_if(path == "", "example dataset not installed")
  spe <- readRDS(path)

  # "hierarchical" is the one that transposes: R holds a block's components as a
  # K_terms x K_groups matrix filled BY ROW from the flat slice, and handing the
  # slice to the shrinkage as it lies has been a bug before.
  for (shrinkage in c("shared", "hierarchical", "half_cauchy")) {
    fit_once <- function() {
      paceModel(spe, celltype_col = "cellType", tau_shrinkage = shrinkage,
                n_iter = 5L, min_iter = 3L, threads = 2L, verbose = FALSE)@fit
    }
    driven <- fit_once()
    reference <- fit_through_r_loop(fit_once())
    expect_identical(loop_outputs(driven), loop_outputs(reference))
  }
})

test_that("the core driver reproduces the R loop on the Gaussian family", {
  # The identity link takes the coefficient-based convergence metric and the
  # residual-variance dispersion step, neither of which the count path reaches.
  set.seed(11)
  n <- 700L
  n_genes <- 16L
  df <- data.frame(celltype = sample(c("A", "B", "C"), n, replace = TRUE),
                   A = runif(n), B = runif(n), C = runif(n),
                   stringsAsFactors = FALSE)
  x_fixed <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  values <- matrix(rnorm(n * n_genes, mean = 4, sd = 1.5), n, n_genes)
  colnames(values) <- paste0("p", seq_len(n_genes))
  counts <- methods::as(methods::as(methods::as(values, "dMatrix"), "generalMatrix"),
                        "CsparseMatrix")
  re_specs <- list(list(group_col = "celltype", formula = ~ 0 + A + B + C))

  fit_once <- function() {
    fit_pace_mvpql_streaming(
      Y = counts, X_fixed = x_fixed, df = df, re_specs = re_specs,
      family = "gaussian", n_iter = 6L, min_iter = 3L, chunk_size = 8L,
      n_threads = 2L, verbose = FALSE)
  }
  driven <- fit_once()
  reference <- fit_through_r_loop(fit_once())
  expect_identical(loop_outputs(driven), loop_outputs(reference))
})

test_that("a binding tau cap warns from the driver exactly as it did from R", {
  skip_if_not_installed("SpatialExperiment")
  path <- system.file("extdata", "bc_xenium_subset.rds", package = "PACE")
  skip_if(path == "", "example dataset not installed")
  spe <- readRDS(path)

  collect_warnings <- function(expression) {
    seen <- character(0)
    withCallingHandlers(force(expression),
                        warning = function(w) {
                          seen <<- c(seen, conditionMessage(w))
                          invokeRestart("muffleWarning")
                        })
    seen
  }
  fit_once <- function() {
    paceModel(spe, celltype_col = "cellType", tau_max = 0.001, n_iter = 4L,
              min_iter = 2L, threads = 2L, verbose = FALSE)@fit
  }
  driven <- NULL
  reference <- NULL
  driven_warnings <- collect_warnings(driven <- fit_once())
  reference_warnings <- collect_warnings(reference <- fit_through_r_loop(fit_once()))

  expect_gt(length(driven_warnings), 0L)
  expect_identical(driven_warnings, reference_warnings)
  expect_identical(loop_outputs(driven), loop_outputs(reference))
})
