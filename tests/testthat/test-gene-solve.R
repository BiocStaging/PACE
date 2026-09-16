# The per-gene penalised weighted least squares solve, against the textbook
# system it is a factored form of. The engine never builds Z; the reference
# below builds it, assembles the whole (p + q) system and solves it with base R,
# so the two implementations share nothing but the maths.

reference_gene_solve <- function(X, terms_list, cell_group_list, blocks, w, z, lam) {
  n <- nrow(X)
  p <- ncol(X)
  q <- sum(vapply(blocks, function(b) b$K_terms * b$K_groups, numeric(1)))
  # Z has one column per (term, group): column col_offset + (t - 1) * K_groups + g
  # carries term t of the cells in group g, and zero elsewhere.
  Z <- matrix(0, n, q)
  for (b in seq_along(blocks)) {
    blk <- blocks[[b]]
    for (t in seq_len(blk$K_terms)) {
      for (g in seq_len(blk$K_groups)) {
        cells <- which(cell_group_list[[b]] == g)
        Z[cells, blk$col_offset + (t - 1L) * blk$K_groups + g] <- terms_list[[b]][cells, t]
      }
    }
  }
  D <- cbind(X, Z)
  out <- list(B = matrix(0, p, ncol(w)), U = matrix(0, q, ncol(w)),
              Ainv_diag = matrix(0, p + q, ncol(w)))
  for (gi in seq_len(ncol(w))) {
    A <- crossprod(D, D * w[, gi])
    A[p + seq_len(q), p + seq_len(q)] <- A[p + seq_len(q), p + seq_len(q)] + diag(lam[, gi], q)
    rhs <- crossprod(D, w[, gi] * z[, gi])
    solution <- solve(A, rhs)
    out$B[, gi] <- solution[seq_len(p)]
    out$U[, gi] <- solution[p + seq_len(q)]
    out$Ainv_diag[, gi] <- diag(solve(A))
  }
  out
}

synthetic_design <- function(n = 600L, n_terms = 4L, n_groups = 5L, n_genes = 6L,
                             second_block = FALSE, seed = 13L) {
  set.seed(seed)
  X <- matrix(1, n, 1L)
  group <- sample.int(n_groups, n, TRUE)
  terms <- cbind(1, matrix(stats::rnorm(n * (n_terms - 1L)), n, n_terms - 1L))
  blocks <- list(list(col_offset = 0L, K_terms = n_terms, K_groups = n_groups))
  terms_list <- list(terms)
  group_list <- list(group)
  if (second_block) {
    image_groups <- 3L
    image <- sample.int(image_groups, n, TRUE)
    blocks[[2]] <- list(col_offset = n_terms * n_groups, K_terms = 2L, K_groups = image_groups)
    terms_list[[2]] <- cbind(1, stats::rnorm(n))
    group_list[[2]] <- image
  }
  q <- sum(vapply(blocks, function(b) b$K_terms * b$K_groups, numeric(1)))
  list(X = X, terms_list = terms_list, group_list = group_list, blocks = blocks,
       cells_by_group = lapply(seq_along(blocks), function(b)
         lapply(seq_len(blocks[[b]]$K_groups), function(g) which(group_list[[b]] == g))),
       w = matrix(abs(stats::rnorm(n * n_genes, 4, 1)), n, n_genes),
       z = matrix(stats::rnorm(n * n_genes), n, n_genes),
       lam = matrix(abs(stats::rnorm(q * n_genes, 2, 0.5)), q, n_genes))
}

solve_with_engine <- function(design, single_precision = FALSE, threads = 1L) {
  PACE:::pace_solve_genes_chunk_cpp(
    design$X, design$w, design$z, design$lam,
    q_total = nrow(design$lam), blocks = design$blocks, terms_list = design$terms_list,
    cells_by_group_list = design$cells_by_group, cell_group_list = design$group_list,
    single_precision = single_precision, n_threads = threads)
}

test_that("the engine solves the same system as an explicit Z", {
  for (second in c(FALSE, TRUE)) {
    design <- synthetic_design(second_block = second)
    fitted <- solve_with_engine(design)
    reference <- reference_gene_solve(design$X, design$terms_list, design$group_list,
                                      design$blocks, design$w, design$z, design$lam)
    label <- if (second) "two blocks" else "one block"
    expect_equal(fitted$B, reference$B, tolerance = 1e-9, label = label)
    expect_equal(fitted$U, reference$U, tolerance = 1e-9, label = label)
    expect_equal(fitted$Ainv_diag, reference$Ainv_diag, tolerance = 1e-9, label = label)
  }
})

test_that("the engine does not depend on the thread count", {
  # The kernel it replaces did: its results moved between one thread and four.
  for (second in c(FALSE, TRUE)) {
    design <- synthetic_design(second_block = second, n_genes = 12L)
    for (single in c(FALSE, TRUE)) {
      expect_identical(solve_with_engine(design, single, 1L),
                       solve_with_engine(design, single, 4L))
    }
  }
})

test_that("a singular gene leaves its column missing rather than a wrong answer", {
  design <- synthetic_design(n_genes = 3L)
  design$w[, 2] <- 0                      # no information at all for this gene
  design$lam[, 2] <- 0                    # and no ridge to identify it
  fitted <- solve_with_engine(design)
  expect_true(all(is.na(fitted$U[, 2])))
  expect_false(any(is.na(fitted$U[, c(1, 3)])))
})

test_that("the engine refuses inconsistent designs", {
  design <- synthetic_design()
  design$group_list[[1]][1] <- 99L
  expect_error(solve_with_engine(design), "index is out of range")
})
