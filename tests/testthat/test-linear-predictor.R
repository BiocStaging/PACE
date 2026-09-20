## The chunk linear predictor and the chunk working response, against the R they
## replaced. These are the two computations that were still in R inside the
## streaming solver's hot loop; nothing else pins them, because the fixture
## gates only see whole fits.

r_eta_block <- function(X_fixed, B, Z, U, genes, p, x1, x1_is_unit) {
  xb <- if (p == 1L) {
    broadcast <- matrix(B[1L, genes], nrow = nrow(X_fixed), ncol = length(genes), byrow = TRUE)
    if (x1_is_unit) broadcast else x1 * broadcast
  } else {
    as.matrix(X_fixed %*% B[, genes, drop = FALSE])
  }
  xb + as.matrix(Z %*% U[, genes, drop = FALSE])
}

test_that("the core's eta block matches the R expression it replaced", {
  skip_if_not_installed("Matrix")
  set.seed(11)
  worst <- c(`1` = 0, `2` = 0, `3` = 0)
  for (n in c(1L, 7L, 300L)) {
    for (q in c(1L, 4L, 23L)) {
      for (p in c(1L, 2L, 3L)) {
        n_genes <- 9L
        Z <- methods::as(Matrix::rsparsematrix(n, q, density = min(1, 3 / max(q, 1)),
                                               rand.x = stats::rnorm), "CsparseMatrix")
        B <- matrix(stats::rnorm(p * n_genes), p, n_genes)
        U <- matrix(stats::rnorm(q * n_genes), q, n_genes)
        ## exact zeros exercise the skipped-column path
        U[sample(length(U), floor(length(U) / 4))] <- 0
        for (x1_is_unit in c(TRUE, FALSE)) {
          x1 <- if (x1_is_unit) rep(1, n) else stats::rnorm(n)
          X_fixed <- if (p == 1L) matrix(x1, n, 1) else
            cbind(x1, matrix(stats::rnorm(n * (p - 1L)), n, p - 1L))
          ## a non-contiguous gene order must work: the chunk is not always 1:m
          for (genes in list(1:3, c(9L, 1L, 5L), 1:9)) {
            reference <- r_eta_block(X_fixed, B, Z, U, genes, p, x1, x1_is_unit)
            for (threads in c(1L, 4L)) {
              got <- pace_eta_block_cpp(if (x1_is_unit) numeric(0) else x1, x1_is_unit,
                                        if (p == 1L) matrix(0, 0, 0) else X_fixed, p,
                                        B, Z, U, as.integer(genes), threads)
              expect_identical(dim(got), dim(reference))
              worst[[as.character(p)]] <- max(worst[[as.character(p)]],
                                              max(abs(got - reference)))
            }
          }
        }
      }
    }
  }
  ## p == 1 is a broadcast multiply on both sides, so it is exact. p > 1 sends
  ## X_fixed %*% B through BLAS in R, whose accumulation is not reproducible
  ## portably; the fixtures are gated at 1e-10, five orders above this.
  expect_identical(worst[["1"]], 0)
  expect_lt(worst[["2"]], 1e-12)
  expect_lt(worst[["3"]], 1e-12)
})

test_that("the core's eta block is invariant to the thread count", {
  skip_if_not_installed("Matrix")
  set.seed(5)
  n <- 400L; q <- 17L; p <- 2L; n_genes <- 12L
  Z <- methods::as(Matrix::rsparsematrix(n, q, density = 0.2, rand.x = stats::rnorm),
                   "CsparseMatrix")
  B <- matrix(stats::rnorm(p * n_genes), p, n_genes)
  U <- matrix(stats::rnorm(q * n_genes), q, n_genes)
  X_fixed <- cbind(1, stats::rnorm(n))
  genes <- as.integer(c(12L, 3L, 7L, 1L))
  one <- pace_eta_block_cpp(numeric(0), FALSE, X_fixed, p, B, Z, U, genes, 1L)
  for (threads in c(2L, 4L, 8L)) {
    expect_identical(pace_eta_block_cpp(numeric(0), FALSE, X_fixed, p, B, Z, U, genes, threads),
                     one)
  }
})

test_that("the chunk working response does not depend on the sub-block width", {
  skip_if_not_installed("Matrix")
  set.seed(7)
  n <- 250L; q <- 6L; p <- 1L; n_genes <- 8L
  counts <- methods::as(Matrix::rsparsematrix(n, n_genes, density = 0.3,
                                              rand.x = function(k) rpois(k, 3) + 1),
                        "CsparseMatrix")
  ambient <- methods::as(Matrix::rsparsematrix(n, n_genes, density = 0.3,
                                               rand.x = function(k) abs(stats::rnorm(k))),
                         "CsparseMatrix")
  Z <- methods::as(Matrix::rsparsematrix(n, q, density = 0.3, rand.x = stats::rnorm),
                   "CsparseMatrix")
  B <- matrix(stats::rnorm(p * n_genes), p, n_genes)
  U <- matrix(stats::rnorm(q * n_genes), q, n_genes)
  genes <- seq_len(n_genes)
  arguments <- list(matrix(0, 0, 0), numeric(0), TRUE, matrix(0, 0, 0), p, B, Z, U,
                    as.integer(genes), counts, ambient, 1L, stats::rnorm(n),
                    runif(n), runif(n_genes, 0.1, 2), numeric(0), FALSE, FALSE, FALSE, n)
  ##                                    sample_weight  nb2  gaussian  seed_iteration
  ## Sub-blocking is a memory strategy: every width must give the same answer.
  reference <- do.call(pace_working_response_chunk_cpp, c(arguments, list(n_genes, 1L)))
  for (sub_genes in c(1L, 2L, 3L, 5L, n_genes)) {
    for (threads in c(1L, 4L)) {
      got <- do.call(pace_working_response_chunk_cpp,
                     c(arguments, list(as.integer(sub_genes), threads)))
      expect_equal(got$z, reference$z, tolerance = 1e-12)
      expect_equal(got$w, reference$w, tolerance = 1e-12)
      expect_equal(got$colsum_w, reference$colsum_w, tolerance = 1e-12)
    }
  }
})


test_that("the Gaussian working response is z = y and a per-gene constant weight", {
  skip_if_not_installed("Matrix")
  set.seed(11)
  n <- 200L; q <- 5L; p <- 1L; n_genes <- 6L
  ## Continuous intensities, not counts: the Gaussian path is for IMC protein.
  dense <- matrix(abs(stats::rnorm(n * n_genes, mean = 2)), n, n_genes)
  counts <- methods::as(methods::as(dense, "Matrix"), "CsparseMatrix")
  ambient <- methods::as(Matrix::rsparsematrix(n, n_genes, density = 0.3,
                                               rand.x = function(k) abs(stats::rnorm(k))),
                         "CsparseMatrix")
  Z <- methods::as(Matrix::rsparsematrix(n, q, density = 0.3, rand.x = stats::rnorm),
                   "CsparseMatrix")
  B <- matrix(stats::rnorm(p * n_genes), p, n_genes)
  U <- matrix(stats::rnorm(q * n_genes), q, n_genes)
  sigma2 <- runif(n_genes, 0.1, 2)
  offset <- stats::rnorm(n)
  rho <- runif(n)

  got <- pace_working_response_chunk_cpp(
    matrix(0, 0, 0), numeric(0), TRUE, matrix(0, 0, 0), p, B, Z, U,
    seq_len(n_genes), counts, ambient, 1L, offset, rho, sigma2, numeric(0),
    FALSE, TRUE, FALSE, n, n_genes, 1L)

  ## z is the raw intensity: it does not depend on the fit, the offset or rho.
  expect_equal(got$z, dense, tolerance = 0)
  ## w is 1/sigma2_g, the same value for every cell of a gene.
  expect_equal(got$w, matrix(rep(1 / sigma2, each = n), n, n_genes), tolerance = 0)
  expect_equal(got$colsum_w, n / sigma2, tolerance = 1e-12)

  ## Neither eta, the offset, rho nor the ambient block may reach the answer.
  other <- pace_working_response_chunk_cpp(
    matrix(0, 0, 0), numeric(0), TRUE, matrix(0, 0, 0), p, B * 3, Z, U * 3,
    seq_len(n_genes), counts, ambient, 1L, offset * 5, rho * 0, sigma2, numeric(0),
    FALSE, TRUE, FALSE, n, n_genes, 1L)
  expect_identical(other$z, got$z)
  expect_identical(other$w, got$w)

  ## Sub-blocking and threading must not move it either.
  for (sub_genes in c(1L, 2L, n_genes)) {
    for (threads in c(1L, 4L)) {
      split <- pace_working_response_chunk_cpp(
        matrix(0, 0, 0), numeric(0), TRUE, matrix(0, 0, 0), p, B, Z, U,
        seq_len(n_genes), counts, ambient, 1L, offset, rho, sigma2, numeric(0),
        FALSE, TRUE, FALSE, n, as.integer(sub_genes), threads)
      expect_identical(split$z, got$z)
      expect_identical(split$w, got$w)
      expect_equal(split$colsum_w, got$colsum_w, tolerance = 1e-12)
    }
  }

  ## A non-positive sigma2 is refused rather than silently giving Inf weights.
  expect_error(
    pace_working_response_chunk_cpp(
      matrix(0, 0, 0), numeric(0), TRUE, matrix(0, 0, 0), p, B, Z, U,
      seq_len(n_genes), counts, ambient, 1L, offset, rho, replace(sigma2, 1, 0),
      numeric(0), FALSE, TRUE, FALSE, n, n_genes, 1L),
    "sigma2")
})


## The blocked eta path. Z is never a general sparse matrix: it is a stack of
## blocks in which every cell belongs to exactly one group, so the same sum can
## be taken a panel of cells at a time with the genes on the inside, reading the
## design once per gene chunk instead of once per gene. That reordering is only
## legitimate if it is EXACT -- the manuscript's fits depend on the bit -- so it
## is gated against the compressed-column walk with identical(), not a tolerance.

## The design of one block, exactly as build_random_design_multi() builds it.
random_design_block <- function(n, n_terms, n_groups, offset) {
  X_terms <- matrix(stats::rnorm(n * n_terms), n, n_terms)
  cell_group <- sample.int(n_groups, n, replace = TRUE)
  list(Z = pace_random_design_block_cpp(X_terms, cell_group - 1L, n_groups)$Z,
       block = list(col_offset = offset, K_terms = n_terms, K_groups = n_groups),
       terms = X_terms, cell_group = cell_group)
}

test_that("the blocked eta path is bit-identical to the compressed-column walk", {
  skip_if_not_installed("Matrix")
  set.seed(19)
  shapes <- list(
    list(c(1L, 1L)),                  # one term, one group: the degenerate block
    list(c(4L, 3L)),                  # one block, as breast cancer has
    list(c(5L, 4L), c(1L, 7L)),       # two blocks, as melanoma has
    list(c(2L, 9L), c(3L, 2L), c(1L, 1L)))
  for (shape in shapes) {
    for (n in c(1L, 13L, 700L)) {
      for (p in c(1L, 2L)) {
        parts <- list(); offset <- 0L
        for (s in shape) {
          parts[[length(parts) + 1L]] <- random_design_block(n, s[1], s[2], offset)
          offset <- offset + s[1] * s[2]
        }
        Z <- do.call(cbind, lapply(parts, `[[`, "Z"))
        q <- offset
        n_genes <- 11L
        B <- matrix(stats::rnorm(p * n_genes), p, n_genes)
        U <- matrix(stats::rnorm(q * n_genes), q, n_genes)
        ## Exact zeros are the one place the two paths could differ: the CSC
        ## walk skips such a column, the blocked path adds 0 * terms instead.
        U[sample(length(U), floor(length(U) / 3))] <- 0
        x1 <- stats::rnorm(n)
        X_fixed <- if (p == 1L) matrix(0, 0, 0) else
          cbind(x1, matrix(stats::rnorm(n * (p - 1L)), n, p - 1L))
        for (genes in list(1:3, c(11L, 2L, 7L, 1L, 4L), 1:11)) {
          for (x1_is_unit in if (p == 1L) c(TRUE, FALSE) else FALSE) {
            reference <- pace_eta_block_cpp(if (x1_is_unit) numeric(0) else x1, x1_is_unit,
                                            X_fixed, p, B, Z, U, as.integer(genes), 1L)
            for (panel in c(0L, 1L, 2L, 64L, 4096L)) {
              for (threads in c(1L, 4L)) {
                got <- pace_eta_block_blocked_cpp(
                  if (x1_is_unit) numeric(0) else x1, x1_is_unit, X_fixed, p, B, Z,
                  lapply(parts, `[[`, "block"), lapply(parts, `[[`, "terms"),
                  lapply(parts, `[[`, "cell_group"), U, as.integer(genes), panel, threads)
                expect_identical(got, reference)
              }
            }
          }
        }
      }
    }
  }
})

test_that("the blocked eta path refuses blocks that are not Z's columns", {
  skip_if_not_installed("Matrix")
  set.seed(23)
  n <- 40L; p <- 1L; n_genes <- 3L
  part <- random_design_block(n, 3L, 4L, 0L)
  Z <- part$Z
  B <- matrix(stats::rnorm(p * n_genes), p, n_genes)
  U <- matrix(stats::rnorm(ncol(Z) * n_genes), ncol(Z), n_genes)
  call_with <- function(block, cell_group = part$cell_group) {
    pace_eta_block_blocked_cpp(numeric(0), TRUE, matrix(0, 0, 0), p, B, Z, list(block),
                               list(part$terms), list(cell_group), U, 1:3, 0L, 1L)
  }
  expect_error(call_with(list(col_offset = 1L, K_terms = 3L, K_groups = 4L)), "tile")
  expect_error(call_with(list(col_offset = 0L, K_terms = 2L, K_groups = 4L)), "n \\* n_terms")
  expect_error(call_with(list(col_offset = 0L, K_terms = 3L, K_groups = 3L)), "cover")
  expect_error(call_with(list(col_offset = 0L, K_terms = 3L, K_groups = 4L),
                         cell_group = replace(part$cell_group, 1L, 9L)),
               "outside its block")
})
