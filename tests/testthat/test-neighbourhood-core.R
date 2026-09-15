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

test_that("kernels match the frNN + sparseMatrix reference exactly", {
  skip_if_not_installed("SpatialExperiment")
  spe <- bc_crop()
  coords <- SpatialExperiment::spatialCoords(spe)
  ct <- spe$cellType
  types <- sort(unique(ct))
  image <- synthetic_images(coords)

  global_new <- PACE:::pace_neighbour_kernel(coords, ct, types, 30, 5, 90, threads = 2L)
  global_ref <- reference_neighbour_kernel(coords, ct, types, 30, 5, 90)
  expect_identical(global_new, global_ref)

  per_image_new <- PACE:::pace_neighbour_kernel(coords, ct, types, 30, 5, 90,
                                                image = image, per_image = TRUE, threads = 2L)
  per_image_ref <- reference_neighbour_kernel(coords, ct, types, 30, 5, 90,
                                              image = image, per_image = TRUE)
  expect_identical(per_image_new, per_image_ref)

  # a type subset: cells of other types are not neighbours but keep their rows
  subset_new <- PACE:::pace_neighbour_kernel(coords, ct, types[1:3], 30, 5, 90)
  subset_ref <- reference_neighbour_kernel(coords, ct, types[1:3], 30, 5, 90)
  expect_identical(subset_new, subset_ref)
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
    expect_identical(new, ref)
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

  expect_identical(buildNeighbourhood(spe, "cellType"),
                   reference_neighbour_kernel(coords, ct, types, 30, 5, 90))

  Y <- t(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  expect_identical(ambientField(spe, "cellType", verbose = FALSE),
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
