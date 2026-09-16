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
                    runif(n), runif(n_genes, 0.1, 2), numeric(0), FALSE, FALSE, n)
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
