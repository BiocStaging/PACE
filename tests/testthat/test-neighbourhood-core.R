# The compiled neighbourhood, ambient-field and anchor code (src/core/) against
# the pure-R implementations it replaces (helper-neighbourhood-reference.R).
# Exact agreement is expected: same distance arithmetic as dbscan::frNN, same
# summation order as Matrix::sparseMatrix and colMeans().

bc_crop <- function(fraction = 0.25) {
  f <- system.file("extdata", "bc_xenium_subset.rds", package = "PACE")
  skip_if(f == "", "example dataset not installed")
  spe <- readRDS(f)
  x <- SpatialExperiment::spatialCoords(spe)[, 1]
  spe[, x <= stats::quantile(x, fraction)]
}

# Two halves of the crop as images, plus one image with a single cell and one
# cell with no image, to exercise the per-image paths.
synthetic_images <- function(coords) {
  image <- ifelse(coords[, 1] <= stats::median(coords[, 1]), "left", "right")
  image[1] <- "single"
  image[2] <- NA_character_
  image
}

## The kernel entry for a cell is a SUM over its neighbours, so it is a
## reduction, and a reduction's result depends on the order the compiler chooses
## to accumulate it. AVX and NEON vectorise it differently: Bioconductor's
## aarch64 builders failed this comparison while every x86_64 build passed, at a
## relative difference of order 1e-16. Exact equality is therefore an assertion
## about the instruction set, not about the kernel.
##
## What still has to match exactly is everything that is not arithmetic -- the
## names, the shapes and which entries are zero, i.e. which cells the search
## found as neighbours at all. A defect in the neighbour search changes the zero
## pattern; a defect in the weighting changes the values by far more than 1e-12.
## Same reasoning as expect_kernel_equal: the ambient field and the exported
## helpers that wrap it are sums over neighbours, so their accumulation order is
## the compiler's choice and differs between AVX and NEON. These return nested
## structures (sparse matrices inside lists), so the comparison is delegated to
## expect_equal(), which walks them and applies the tolerance to the numeric
## leaves while still requiring the shapes, names and classes to match.
expect_numeric_equal <- function(new, ref, tol = 1e-12) {
  expect_equal(new, ref, tolerance = tol)
}

expect_kernel_equal <- function(new, ref, tol = 1e-12) {
  expect_identical(names(new), names(ref))
  for (nm in names(new)) {
    expect_identical(dim(new[[nm]]), dim(ref[[nm]]))
    expect_identical(dimnames(new[[nm]]), dimnames(ref[[nm]]))
    a <- as.numeric(new[[nm]]); b <- as.numeric(ref[[nm]])
    expect_identical(a == 0, b == 0)
    rel <- abs(a - b) / pmax(abs(b), 1e-300)
    expect_lt(max(rel[is.finite(rel)], 0), tol)
  }
}

test_that("kernels match the frNN + sparseMatrix reference", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop()
  coords <- SpatialExperiment::spatialCoords(spe)
  ct <- spe$cellType
  types <- sort(unique(ct))
  image <- synthetic_images(coords)

  global_new <- PACE:::pace_neighbour_kernel(coords, ct, types, 30, 5, 90, threads = 2L)
  global_ref <- reference_neighbour_kernel(coords, ct, types, 30, 5, 90)
  expect_kernel_equal(global_new, global_ref)

  per_image_new <- PACE:::pace_neighbour_kernel(coords, ct, types, 30, 5, 90,
                                                image = image, per_image = TRUE, threads = 2L)
  per_image_ref <- reference_neighbour_kernel(coords, ct, types, 30, 5, 90,
                                              image = image, per_image = TRUE)
  expect_kernel_equal(per_image_new, per_image_ref)

  # a type subset: cells of other types are not neighbours but keep their rows
  subset_new <- PACE:::pace_neighbour_kernel(coords, ct, types[1:3], 30, 5, 90)
  subset_ref <- reference_neighbour_kernel(coords, ct, types[1:3], 30, 5, 90)
  expect_kernel_equal(subset_new, subset_ref)
})

test_that("neighbour search reproduces frNN at exactly eps, with duplicates and tiny images", {
  skip_if_not_installed("dbscan")
  set.seed(11)
  frnn_counts <- function(coords, eps) lengths(dbscan::frNN(coords, eps = eps)$id)

  # integer lattice: many pairs at exactly eps, and bin edges on lattice lines
  lattice <- as.matrix(expand.grid(x = 0:40, y = 0:40))
  for (eps in c(1, 2, 5)) {
    expect_identical(PACE:::pace_neighbour_counts_cpp(lattice + 0, rep(-1L, nrow(lattice)),
                                                      FALSE, eps, 2L),
                     frnn_counts(lattice + 0, eps))
  }

  # pairs placed at distance eps along 3-4-5 directions from random origins, so
  # the rounded squared distance falls on, just under or just over eps^2
  n_pairs <- 25000
  eps <- 7.3
  origin <- cbind(stats::runif(n_pairs, -1e3, 1e3), stats::runif(n_pairs, -1e3, 1e3))
  direction <- cbind(sample(c(-0.6, 0.6, -0.8, 0.8), n_pairs, TRUE), 0)
  direction[, 2] <- ifelse(abs(direction[, 1]) == 0.6, 0.8, 0.6) * sample(c(-1, 1), n_pairs, TRUE)
  partner <- origin + eps * direction
  adversarial <- rbind(origin, partner)
  expect_identical(PACE:::pace_neighbour_counts_cpp(adversarial, rep(-1L, nrow(adversarial)),
                                                    FALSE, eps, 4L),
                   frnn_counts(adversarial, eps))

  # points on the grid origin offsets: x = xmin + k * eps exactly and just beside it
  edge_x <- rep(c(0, 7.3, 14.6, 21.9), each = 3) + c(0, 1e-12, -1e-12)
  edge <- cbind(c(edge_x, edge_x), c(rep(0, 12), rep(7.3, 12)))
  expect_identical(PACE:::pace_neighbour_counts_cpp(edge, rep(-1L, nrow(edge)), FALSE, 7.3, 1L),
                   frnn_counts(edge, 7.3))

  # duplicates at distance 0 are neighbours of each other
  dup <- rbind(c(1, 1), c(1, 1), c(1, 1), c(3, 1))
  expect_identical(PACE:::pace_neighbour_counts_cpp(dup, rep(-1L, 4), FALSE, 2, 1L),
                   frnn_counts(dup, 2))

  # per image: images with fewer than 2 cells, and cells with no image, get none
  group <- c(0L, 0L, 1L, -1L)
  expect_identical(PACE:::pace_neighbour_counts_cpp(dup, group, TRUE, 2, 1L), c(1L, 1L, 0L, 0L))
  expect_identical(PACE:::pace_neighbour_counts_cpp(matrix(numeric(0), 0, 2), integer(0), FALSE, 2, 2L),
                   integer(0))
})

test_that("edge correction matches the quadrature reference, interior cells included", {
  skip_if_not_installed("SpatialExperiment")
  coords <- SpatialExperiment::spatialCoords(bc_crop())
  expect_identical(PACE:::pace_area_fraction(coords, 15), reference_area_fraction(coords, 15))
  expect_identical(PACE:::pace_area_fraction(coords, 90), reference_area_fraction(coords, 90))
  corners <- rbind(c(0, 0), c(10, 0), c(0, 10), c(10, 10), c(5, 5), c(5, 0))
  expect_identical(PACE:::pace_area_fraction(corners, 3), reference_area_fraction(corners, 3))
  single <- matrix(c(2, 2), 1, 2)
  expect_identical(PACE:::pace_area_fraction(single, 3), reference_area_fraction(single, 3))
})

test_that("ambient field W matches the reference, per image and without edge correction", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop()
  coords <- SpatialExperiment::spatialCoords(spe)
  Y <- t(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  ct <- spe$cellType
  types <- sort(unique(ct))
  image <- factor(synthetic_images(coords))
  for (edge_correct in c(TRUE, FALSE)) {
    new <- PACE:::pace_ambient_field(coords, PACE:::.pace_as_dgc(Y), ct, image, types, 5,
                                     edge_correct = edge_correct, verbose = FALSE, threads = 2L)
    ref <- reference_ambient_field(coords, Y, ct, image, types, 5, edge_correct = edge_correct)
    expect_numeric_equal(new, ref)
  }
})

test_that("anchors and detection rates match the reference", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop()
  coords <- SpatialExperiment::spatialCoords(spe)
  Y <- t(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  ct <- spe$cellType
  types <- sort(unique(ct))
  image <- synthetic_images(coords)
  image[3:10] <- "small"   # an image below the 50-cell threshold
  new <- PACE:::pace_anchors(coords, PACE:::.pace_as_dgc(Y), factor(ct, levels = types), image,
                             types, verbose = FALSE, threads = 2L)
  ref <- reference_anchors(coords, Y, factor(ct, levels = types), image, types)
  expect_identical(new$mask, ref$mask)
  expect_identical(new$idx, ref$idx)

  focals <- c(types, "Absent_Type")
  new_rate <- PACE:::pace_group_column_means_cpp(PACE:::.pace_as_dgc(Y),
                                                 PACE:::.pace_codes(ct, focals), length(focals),
                                                 detection = TRUE, n_threads = 2L)
  ref_rate <- reference_detection_rate(Y, ct, focals)
  expect_identical(new_rate[seq_along(types), ], unname(ref_rate[seq_along(types), ]))
  expect_true(all(is.nan(new_rate[length(focals), ])))
})

test_that("results do not depend on the number of threads", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop(0.5)
  coords <- SpatialExperiment::spatialCoords(spe)
  Y <- PACE:::.pace_as_dgc(t(as.matrix(SummarizedExperiment::assay(spe, "counts"))))
  ct <- spe$cellType
  types <- sort(unique(ct))
  image <- synthetic_images(coords)
  kernel <- function(threads, per_image) {
    PACE:::pace_neighbour_kernel(coords, ct, types, 30, 5, 90, image = image,
                                 per_image = per_image, threads = threads)
  }
  expect_identical(kernel(1L, FALSE), kernel(4L, FALSE))
  expect_identical(kernel(1L, TRUE), kernel(4L, TRUE))
  ambient <- function(threads) {
    PACE:::pace_ambient_field(coords, Y, ct, factor(image), types, 5, verbose = FALSE,
                              validate = FALSE, threads = threads)
  }
  expect_identical(ambient(1L), ambient(4L))
  expect_identical(PACE:::pace_area_fraction(coords, 15, threads = 1L),
                   PACE:::pace_area_fraction(coords, 15, threads = 4L))
  anchors <- function(threads) {
    PACE:::pace_anchors(coords, Y, factor(ct, levels = types), image, types,
                        verbose = FALSE, threads = threads)
  }
  expect_identical(anchors(1L), anchors(4L))
})

test_that("exported neighbourhood helpers keep their outputs", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop()
  coords <- SpatialExperiment::spatialCoords(spe)
  ct <- as.character(spe$cellType)
  types <- sort(unique(ct))

  expect_numeric_equal(buildNeighbourhood(spe, "cellType"),
                   reference_neighbour_kernel(coords, ct, types, 30, 5, 90))

  Y <- t(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  expect_numeric_equal(ambientField(spe, "cellType", verbose = FALSE),
                   reference_ambient_field(coords, Y, ct, factor(rep("all", nrow(Y))), types, 5))

  fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
  full <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
  expect_s3_class(suppressMessages(suppressWarnings(anchorGenes(fit, full))), "data.frame")
  full_ct <- as.character(full$cellType)
  full_types <- sort(unique(full_ct))
  full_Y <- t(as.matrix(SummarizedExperiment::assay(full, "counts")))
  ref_mask <- suppressWarnings(reference_anchors(SpatialExperiment::spatialCoords(full), full_Y,
                                                 full_ct, rep("all", nrow(full_Y)), full_types)$mask)
  ref_means <- vapply(full_types, function(tt) colMeans(full_Y[full_ct == tt, , drop = FALSE]),
                      numeric(ncol(full_Y)))
  ref_owner <- apply(ref_means, 1, max)
  names(ref_owner) <- colnames(full_Y)
  ref_rows <- lapply(full_types, function(focal) {
    ag <- colnames(full_Y)[ref_mask[focal, ] == 1]
    ranked <- ag[order(ref_owner[ag], decreasing = TRUE)]
    data.frame(cellType = focal, n_anchors = length(ag),
               top_anchors = paste(utils::head(ranked, 10L), collapse = ", "),
               stringsAsFactors = FALSE)
  })
  expect_identical(suppressMessages(suppressWarnings(anchorGenes(fit, full))),
                   do.call(rbind, ref_rows))
})

test_that("cells without a modelled cell type are refused with a clear error", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop()
  types <- setdiff(sort(unique(spe$cellType)), "Dendritic_Cell")
  expect_error(paceModel(spe, celltype_col = "cellType", types = types, n_iter = 1L,
                         threads = 1L, verbose = FALSE),
               "not in `types`")
})

test_that("same-type fractions match the reference, including cells with no neighbours", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop()
  coords <- SpatialExperiment::spatialCoords(spe)
  Y <- t(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  ct <- spe$cellType
  types <- sort(unique(ct))
  image <- synthetic_images(coords)
  image[3:10] <- "small"
  coords[11, ] <- c(1e6, 1e6)   # isolated cell inside a large image: fraction 1
  labels <- as.character(ct)
  new_fraction <- PACE:::pace_same_type_fraction_cpp(coords, PACE:::.pace_codes(labels, unique(labels)),
                                                     PACE:::.pace_image_codes(image)$code,
                                                     radius = 30, min_image_cells = 50L, n_threads = 2L)
  ref_fraction <- reference_anchors(coords, Y, factor(ct, levels = types), image, types)$same_frac
  expect_identical(new_fraction, ref_fraction)
  expect_identical(new_fraction[11], 1)
  expect_true(all(is.na(new_fraction[3:10])))
})

test_that("neighbour sets and distances equal frNN on boundary layouts", {
  skip_if_not_installed("dbscan")
  set.seed(12)
  ## Returns the points where the two neighbour sets disagree in a way that is
  ## NOT explained by the radius boundary; character(0) means they agree.
  ##
  ## Exact set equality is the wrong assertion at the boundary. A pair whose
  ## separation is indistinguishable from eps in double precision can fall
  ## either side of the radius depending on how each package's C++ was
  ## compiled: PACE suppresses FMA contraction (src/core/fp_no_contract.hpp)
  ## but dbscan's build is outside this package's control and differs by
  ## platform. The `pairs` layout below is built from 0.6-0.8-1 triangles
  ## scaled by eps, so 1,758 of its 3,000 separations land exactly ON the
  ## radius and 2,207 within two ulps of it; requiring identical sets there
  ## tests the two toolchains against each other, not the neighbour search.
  ## Bioconductor's builders failed all six platforms on precisely that.
  ##
  ## So a neighbour found by only one side is accepted only when its distance
  ## sits within a few ulps of eps. Everything else must still match: any
  ## disagreement away from the boundary is a real defect, and the neighbours
  ## both sides found must agree on distance to the same few ulps.
  compare_sets <- function(coords, eps) {
    lists <- PACE:::pace_neighbour_lists_cpp(coords, rep(-1L, nrow(coords)), FALSE, eps, 2L)
    reference <- dbscan::frNN(coords, eps = eps)
    offsets <- lists$offsets
    slack <- 8 * .Machine$double.eps * eps
    on_radius <- function(i, j) {
      d <- sqrt(sum((coords[i, ] - coords[j, ])^2))
      abs(d - eps) <= slack
    }
    problems <- character(0)
    for (i in seq_len(nrow(coords))) {
      span <- if (offsets[i + 1] > offsets[i]) (offsets[i] + 1):offsets[i + 1] else integer(0)
      ours <- lists$neighbours[span]
      ours_dist <- lists$distances[span]
      ref_order <- order(reference$id[[i]])
      theirs <- reference$id[[i]][ref_order]
      theirs_dist <- reference$dist[[i]][ref_order]
      if (identical(ours, theirs) && identical(ours_dist, theirs_dist)) next

      disputed <- c(setdiff(ours, theirs), setdiff(theirs, ours))
      off_boundary <- disputed[!vapply(disputed, function(j) on_radius(i, j), logical(1))]
      if (length(off_boundary)) {
        problems <- c(problems, sprintf(
          "point %d: neighbour(s) %s disagree away from the radius",
          i, paste(off_boundary, collapse = ", ")))
        next
      }
      shared <- intersect(ours, theirs)
      if (length(shared)) {
        gap <- abs(ours_dist[match(shared, ours)] - theirs_dist[match(shared, theirs)])
        if (any(gap > slack)) {
          problems <- c(problems, sprintf(
            "point %d: shared neighbour distances differ by %.3g", i, max(gap)))
        }
      }
    }
    problems
  }
  lattice <- as.matrix(expand.grid(0:25, 0:25)) + 0
  expect_equal(compare_sets(lattice, 1), character(0))
  expect_equal(compare_sets(lattice, 2), character(0))
  origin <- cbind(stats::runif(3000, 0, 300), stats::runif(3000, 0, 300))
  pairs <- rbind(origin, origin + 7.3 * cbind(0.6, 0.8)[rep(1, 3000), ])
  expect_equal(compare_sets(pairs, 7.3), character(0))
  edge_x <- rep(c(0, 15, 30, 45), each = 3) * c(1, 1 + .Machine$double.eps, 1 - .Machine$double.eps)
  edge <- cbind(c(edge_x, edge_x), c(rep(0, 12), rep(15, 12)))
  expect_equal(compare_sets(edge, 15), character(0))
  duplicates <- rbind(c(2, 2), c(2, 2), c(2, 3), c(9, 9))
  expect_equal(compare_sets(duplicates, 1), character(0))
})

# Calls that hung or over-allocated before input validation (zero, denormal or
# huge radii) run in a separate Rscript process with a timeout, so a regression
# fails the test instead of hanging the whole suite. Base system2() is used, not
# callr: processx's SIGCHLD handler would reap the forked workers of later fits.
run_isolated <- function(fun, timeout = 60) {
  script <- tempfile(fileext = ".R")
  result_file <- tempfile(fileext = ".rds")
  on.exit(unlink(c(script, result_file)), add = TRUE)
  writeLines(c(paste0("isolated <- ", paste(deparse(fun), collapse = "\n")),
               sprintf("saveRDS(isolated(), %s)", deparse(result_file))), script)
  status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                                     stdout = FALSE, stderr = FALSE, timeout = timeout,
                                     env = paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = ":")))))
  if (!file.exists(result_file))
    stop("isolated R process did not finish (exit status ", status, "; 124 means timeout)")
  readRDS(result_file)
}

test_that("invalid radii and bandwidths are refused quickly instead of hanging", {
  skip_if_not_installed("SpatialExperiment")
  messages <- run_isolated(function() {
    spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
    x <- SpatialExperiment::spatialCoords(spe)[, 1]
    spe <- spe[, x <= stats::quantile(x, 0.1)]
    coords <- SpatialExperiment::spatialCoords(spe)
    ct <- spe$cellType
    types <- sort(unique(ct))
    Y <- PACE:::.pace_as_dgc(t(as.matrix(SummarizedExperiment::assay(spe, "counts"))))
    image <- factor(rep("all", nrow(coords)))
    n <- nrow(coords)
    caught <- function(expr) tryCatch({ force(expr); "NO ERROR" }, error = function(e) conditionMessage(e))
    out <- list()
    for (bad in list(0, -5, NaN, Inf, NA_real_)) {
      key <- format(bad)
      out[[paste("kernel h_tech", key)]] <- caught(PACE:::pace_neighbour_kernel(coords, ct, types, 30, bad, 90))
      out[[paste("kernel h_bio", key)]] <- caught(PACE:::pace_neighbour_kernel(coords, ct, types, bad, 5, 90))
      out[[paste("kernel eps", key)]] <- caught(PACE:::pace_neighbour_kernel(coords, ct, types, 30, 5, bad))
      out[[paste("ambient h_tech", key)]] <- caught(PACE:::pace_ambient_field(coords, Y, ct, image, types, bad,
                                                                             verbose = FALSE))
      out[[paste("area r", key)]] <- caught(PACE:::pace_area_fraction(coords, bad))
      out[[paste("core counts eps", key)]] <- caught(PACE:::pace_neighbour_counts_cpp(coords + 0, rep(-1L, n),
                                                                                     FALSE, bad, 1L))
      out[[paste("core kernels h_tech", key)]] <- caught(PACE:::pace_neighbour_kernels_cpp(
        coords + 0, rep(0L, n), 1L, rep(-1L, n), FALSE, 30, bad, 90, 1L))
      out[[paste("core ambient h_tech", key)]] <- caught(PACE:::pace_ambient_field_cpp(
        coords + 0, rep(0L, n), rep(0L, n), 1L, bad, TRUE, 1, 0, 1L))
    }
    out[["paceModel h_tech 0"]] <- caught(PACE::paceModel(spe, celltype_col = "cellType", h_tech = 0,
                                                         n_iter = 1L, threads = 1L, verbose = FALSE))
    out[["ambientField h_tech 0"]] <- caught(PACE::ambientField(spe, "cellType", h_tech = 0, verbose = FALSE))
    # a finite h_tech whose ambient radius 3 * h_tech overflows
    out[["ambient h_tech huge"]] <- caught(PACE:::pace_ambient_field(coords, Y, ct, image, types, 1e308,
                                                                    verbose = FALSE))
    out[["core ambient h_tech huge"]] <- caught(PACE:::pace_ambient_field_cpp(
      coords + 0, rep(0L, n), rep(0L, n), 1L, 1e308, TRUE, 1, 0, 1L))
    out[["paceModel h_tech huge"]] <- caught(PACE::paceModel(spe, celltype_col = "cellType", h_tech = 1e308,
                                                            n_iter = 1L, threads = 1L, verbose = FALSE))
    out
  })
  for (name in names(messages)) {
    parameter <- sub("^.* (h_tech|h_bio|eps|r) .*$", "\\1", name)
    expect_false(identical(messages[[name]], "NO ERROR"), info = name)
    expect_match(messages[[name]], if (grepl("huge", name)) "3 \\* h_tech" else parameter, info = name)
  }
})

test_that("a denormal radius over a huge extent stays exact without a huge grid", {
  skip_if_not_installed("dbscan")
  result <- run_isolated(function() {
    coords <- rbind(c(0, 0), c(0, 0), c(1e20, 3e20), c(1e20, 3e20), c(5, 5), c(-3e20, 2))
    list(counts = PACE:::pace_neighbour_counts_cpp(coords, rep(-1L, nrow(coords)), FALSE, 1e-320, 2L),
         reference = lengths(dbscan::frNN(coords, eps = 1e-320)$id),
         small_extent = PACE:::pace_neighbour_counts_cpp(coords[c(1, 2, 5), ], rep(-1L, 3), FALSE, 1e-320, 1L))
  })
  expect_identical(result$counts, result$reference)
  expect_identical(result$small_extent, c(1L, 1L, 0L))
})

test_that("non-finite coordinates, a third coordinate column and non-finite counts are refused", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop(0.1)
  coords <- SpatialExperiment::spatialCoords(spe)
  ct <- spe$cellType
  types <- sort(unique(ct))
  for (bad_value in c(Inf, -Inf, NaN, NA)) {
    bad <- coords
    bad[5, 1] <- bad_value
    expect_error(PACE:::pace_neighbour_kernel(bad, ct, types, 30, 5, 90), "finite")
    expect_error(PACE:::pace_neighbour_counts_cpp(bad + 0, rep(-1L, nrow(bad)), FALSE, 90, 1L),
                 "finite")
  }
  expect_error(PACE:::pace_neighbour_kernel(cbind(coords, z = 0), ct, types, 30, 5, 90),
               "exactly two columns")
  Y <- t(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  Y[1, 1] <- NA
  expect_error(PACE:::pace_anchors(coords, Y, ct, rep("all", nrow(Y)), types, verbose = FALSE),
               "finite")
  counts <- PACE:::.pace_as_dgc(t(as.matrix(SummarizedExperiment::assay(spe, "counts"))))
  counts@x[1] <- NaN
  expect_error(PACE:::pace_group_column_means_cpp(counts, rep(0L, nrow(counts)), 1L, TRUE, 1L), "finite")
})

test_that("DelayedArray assays work in anchorGenes() and ambientField()", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("DelayedArray")
  spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
  delayed <- spe
  SummarizedExperiment::assay(delayed, "counts", withDimnames = FALSE) <-
    DelayedArray::DelayedArray(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  expect_s4_class(SummarizedExperiment::assay(delayed, "counts"), "DelayedMatrix")
  fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
  expect_identical(suppressMessages(suppressWarnings(anchorGenes(fit, delayed))),
                   suppressMessages(suppressWarnings(anchorGenes(fit, spe))))
  expect_identical(ambientField(delayed, "cellType", verbose = FALSE),
                   ambientField(spe, "cellType", verbose = FALSE))
})

test_that("paceModel refuses spatial coordinates with a third column", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop(0.1)
  coords <- SpatialExperiment::spatialCoords(spe)
  spe3 <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = SummarizedExperiment::assay(spe, "counts")),
    colData = SummarizedExperiment::colData(spe),
    spatialCoords = cbind(coords, z = 1))
  expect_error(paceModel(spe3, celltype_col = "cellType", n_iter = 1L, threads = 1L, verbose = FALSE),
               "exactly two columns")
})

test_that("the W identity check passes when no cell has a heterotypic neighbour", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop(0.25)
  coords <- SpatialExperiment::spatialCoords(spe)
  Y <- t(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  ct <- as.character(spe$cellType)
  types <- sort(unique(ct))
  n <- nrow(Y)
  image_all <- factor(rep("all", n))
  # The expected field is an n x n W with no entries. (The pure-R reference cannot
  # build it: sparseMatrix() fails on the empty triplet list.)
  empty_field <- function(image) {
    list(W = Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n)),
         image_idx = as.integer(image), n_images = nlevels(image))
  }

  # a technical radius so small that no two cells are within 3 * h_tech
  tiny <- ambientField(spe, "cellType", h_tech = 0.01, verbose = FALSE)
  expect_identical(tiny, empty_field(image_all))
  expect_identical(ambientField(spe, "cellType", h_tech = 0.1, verbose = FALSE), empty_field(image_all))

  # a single cell type: every neighbour is homotypic
  one_type <- rep("Tumour", n)
  single <- PACE:::pace_ambient_field(coords, PACE:::.pace_as_dgc(Y), one_type, image_all, "Tumour", 5,
                                      verbose = FALSE)
  expect_identical(single, empty_field(image_all))

  # one cell per image: no image has two cells
  per_cell <- factor(seq_len(n))
  isolated <- PACE:::pace_ambient_field(coords, PACE:::.pace_as_dgc(Y), ct, per_cell, types, 5,
                                        verbose = FALSE)
  expect_identical(isolated, empty_field(per_cell))
})
