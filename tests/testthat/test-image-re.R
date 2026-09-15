# Image-level neighbour slopes on a kernel column that is constant over all
# cells. With only 33 B cells in this crop, the sparse-pair drop zeroes the
# B_Cell column for every focal type; scale() then returned NaN for every cell
# and model.matrix() dropped all rows, so the fit stopped.

test_that("image slopes leave out a kernel column that is constant over all cells", {
  skip_if_not_installed("SpatialExperiment")
  f <- system.file("extdata", "mel_cosmx_subset.rds", package = "PACE")
  skip_if(f == "", "example dataset not installed")
  spe <- readRDS(f)
  images <- c("32158_19", "32176_5", "32177_4", "32178_7", "32155_16", "32181_11", "32180_9")
  spe <- spe[, spe$image %in% images]

  for (option in c("slopes", "condition_slopes")) {
    expect_message(
      fit <- paceModel(spe, celltype_col = "cellType", condition_col = "Responder",
                       image_col = "image", kernel_per_image = TRUE,
                       image_re = option, n_iter = 2L, threads = 2L, verbose = TRUE),
      "no image slope for constant kernel column\\(s\\) B_Cell")
    image_terms <- fit@fit$re_meta$blocks[[2]]$term_levels
    expect_false(any(grepl("B_Cell", image_terms)))
    expect_true("Tumour_imgz" %in% image_terms)
    expect_true(all(is.finite(fit@fit$U)))
  }

  # every column constant: the image block keeps only its intercept
  for (option in c("slopes", "condition_slopes")) {
    fit <- suppressMessages(
      paceModel(spe, celltype_col = "cellType", condition_col = "Responder",
                image_col = "image", kernel_per_image = TRUE, drop_sparse_neff = 1e9,
                image_re = option, n_iter = 2L, threads = 2L, verbose = TRUE))
    expect_identical(fit@fit$re_meta$blocks[[2]]$term_levels, "(Intercept)")
  }
})

test_that("condition image slopes without a condition are refused", {
  spe <- readRDS(system.file("extdata", "mel_cosmx_subset.rds", package = "PACE"))
  expect_error(
    paceModel(spe[, 1:500], celltype_col = "cellType", image_col = "image",
              image_re = "condition_slopes", verbose = FALSE),
    "needs a `condition_col`")
})
