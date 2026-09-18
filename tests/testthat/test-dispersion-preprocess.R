# The optimisers and the preprocessing that used to be R calls, against the R
# they replace. The special functions themselves are still R's, passed into the
# core, so these comparisons are expected to be exact.

test_that("the dispersion MLE reproduces optimize() over dnbinom()", {
  set.seed(4)
  n <- 3000L
  mu <- exp(stats::rnorm(n, -1, 1))
  y <- stats::rnbinom(n, size = mu / 0.7, mu = mu)
  nb1 <- function(log_alpha) {
    a <- exp(log_alpha)
    -sum(stats::dnbinom(y, size = mu / a, mu = mu, log = TRUE))
  }
  nb2 <- function(log_alpha) {
    a <- exp(log_alpha)
    -sum(stats::dnbinom(y, size = 1 / a, mu = mu, log = TRUE))
  }
  expect_identical(PACE:::pace_dispersion_mle_cpp(y, mu, FALSE, FALSE, Inf, FALSE),
                   exp(stats::optimize(nb1, interval = c(-6, 4))$minimum))
  expect_identical(PACE:::pace_dispersion_mle_cpp(y, mu, TRUE, FALSE, Inf, FALSE),
                   exp(stats::optimize(nb2, interval = c(-6, 4))$minimum))
  # too few cells, and counts that are not whole numbers, have no MLE
  expect_true(is.na(PACE:::pace_dispersion_mle_cpp(y[1:5], mu[1:5], FALSE, FALSE, Inf, FALSE)))
  # the zero-count collapse is the same algebra, so it lands within the
  # optimiser's own tolerance of the full sum
  collapsed <- PACE:::pace_dispersion_mle_cpp(y, mu, FALSE, TRUE, Inf, FALSE)
  expect_equal(collapsed, PACE:::pace_dispersion_mle_cpp(y, mu, FALSE, FALSE, Inf, FALSE),
               tolerance = 1e-4)
})

test_that("the closed-form NB1 objective finds the same alpha as the density it replaces", {
  # alpha_fast_density minimises the NB1 log-likelihood written out, instead of
  # calling Rf_dnbinom_mu once per cell per Brent step. Once size = mu/a the two
  # per-cell logs become scalars, lgamma(x + 1) drops out (it cannot move the
  # argmin), and small integer counts use sum log(s + k) rather than differencing
  # two lgammas. It is 5x faster and is the default, so it needs to agree with
  # the density it replaced across the range of alpha the fits actually reach.
  set.seed(11)
  n <- 20000L
  for (alpha_true in c(0.05, 0.3, 1, 3, 20)) {
    mu <- stats::rgamma(n, 2, 2) * stats::runif(n, 0.2, 5)
    y <- stats::rnbinom(n, size = mu / alpha_true, mu = mu)
    slow <- PACE:::pace_dispersion_mle_cpp(y, mu, FALSE, TRUE, Inf, FALSE)
    fast <- PACE:::pace_dispersion_mle_cpp(y, mu, FALSE, TRUE, Inf, TRUE)
    # Brent's own tolerance on log alpha is DBL_EPSILON^0.25, about 1.2e-4, so
    # anything at 1e-6 is far inside the optimiser's resolution.
    expect_equal(fast, slow, tolerance = 1e-6,
                 label = sprintf("alpha_true %g", alpha_true))
  }
  # Counts above the small-count limit take the lgamma branch rather than the
  # sum-of-logs one, so give it a gene whose counts are mostly well above it.
  # (Non-integer counts would exercise the same branch, but R's dnbinom warns on
  # every one of them, and the density is called once per cell per Brent step --
  # dispersion_chunk counts such genes separately for exactly that reason.)
  mu <- stats::rgamma(n, 2, 2) * 200
  y <- stats::rnbinom(n, size = mu / 0.4, mu = mu)
  expect_gt(stats::median(y), 8)
  expect_equal(PACE:::pace_dispersion_mle_cpp(y, mu, FALSE, TRUE, Inf, TRUE),
               PACE:::pace_dispersion_mle_cpp(y, mu, FALSE, TRUE, Inf, FALSE),
               tolerance = 1e-6)
  # nb2 has no closed form here and must ignore the flag entirely
  y <- stats::rnbinom(n, size = 1 / 0.5, mu = mu)
  expect_identical(PACE:::pace_dispersion_mle_cpp(y, mu, TRUE, FALSE, Inf, TRUE),
                   PACE:::pace_dispersion_mle_cpp(y, mu, TRUE, FALSE, Inf, FALSE))
})

test_that("the prior degrees of freedom reproduce uniroot() over trigamma()", {
  reference <- function(v_obs, n_used, d0_max = 5) {
    if (n_used < 5L || !is.finite(v_obs)) return(d0_max)
    excess <- v_obs - trigamma(0.5)
    if (excess <= 0) return(d0_max)
    root <- tryCatch(stats::uniroot(function(x) trigamma(x) - excess,
                                    lower = 1e-4, upper = 1e4)$root,
                     error = function(e) NA_real_)
    if (!is.finite(root)) return(d0_max)
    min(2 * root, d0_max)
  }
  for (v in c(0.1, 1.5, 4.9, 5.2, 12, 1e3, 1e-6, NA_real_)) {
    expect_identical(PACE:::pace_estimate_d0_cpp(v, 900, 5), reference(v, 900))
  }
  expect_identical(PACE:::pace_estimate_d0_cpp(2, 4, 5), 5)
})

test_that("the kernel preprocessing reproduces the R loops it replaces", {
  set.seed(9)
  n <- 400L
  types <- c("a", "b", "c")
  celltype <- sample(types, n, TRUE)
  image <- sample(c("i1", "i2"), n, TRUE)
  kernel <- matrix(abs(stats::rnorm(n * 3)), n, 3, dimnames = list(NULL, types))
  kernel[celltype == "b", 3] <- 1                    # constant for one focal

  reference_drop <- function(K, neff_min) {
    for (focal in types) {
      cells <- which(celltype == focal)
      for (nb in seq_along(types)) {
        centred <- K[cells, nb] - mean(K[cells, nb], na.rm = TRUE)
        max_sq <- max(centred^2, na.rm = TRUE)
        n_eff <- if (max_sq > 0) sum(centred^2, na.rm = TRUE) / max_sq else 0
        if (n_eff < neff_min) K[cells, nb] <- 0
      }
    }
    K
  }
  expect_identical(PACE:::pace_drop_sparse_k(kernel, celltype, types, 30, verbose = FALSE),
                   reference_drop(kernel, 30))

  reference_centre <- function(K) {
    for (t in seq_along(types)) K[, t] <- K[, t] - stats::ave(K[, t], image, celltype, FUN = mean)
    K
  }
  expect_identical(PACE:::pace_center_within_image(kernel, celltype, image, types),
                   reference_centre(kernel))

  standardised <- PACE:::pace_standardise_cpp(kernel[, 1])
  expect_identical(standardised$values, as.numeric(scale(kernel[, 1])))
  # the divisor is scale()'s own -- the root mean square of the centred column,
  # which differs from stats::sd() in the last bit because sd() re-centres with
  # the refined mean. Both are zero exactly when the column is constant, which
  # is all the caller asks of it.
  expect_equal(standardised$sd, stats::sd(kernel[, 1]), tolerance = 1e-14)
  constant <- PACE:::pace_standardise_cpp(rep(2, 50))
  expect_identical(constant$sd, 0)
})

test_that("the anchor decision reproduces the R comparison it replaces", {
  set.seed(15)
  types <- c("a", "b", "c", "d")
  core_means <- matrix(abs(stats::rnorm(4 * 30, 0.2, 0.3)), 4, 30, dimnames = list(types, NULL))
  owner_mean <- apply(core_means, 2, max)
  owner <- types[apply(core_means, 2, which.max)]
  reference <- t(vapply(types, function(x)
    as.numeric(owner != x & owner_mean > 0.1 & (core_means[x, ] / pmax(owner_mean, 1e-9) < 0.1)),
    numeric(30)))
  fitted <- PACE:::pace_anchor_mask_cpp(core_means, 0.1, 0.1)
  expect_identical(unname(fitted$mask), unname(reference))
  expect_identical(fitted$n_anchor, as.integer(rowSums(reference)))
})
