# plotResponseCurve() decided which arm carried the responder interaction from
# levels(d$arm)[1], the ALPHABETICALLY first level, while the model's reference
# is whatever model.matrix() used -- which the fit records in params$resp_term.
# The two agree only by coincidence, and when they disagree BOTH dashed slopes
# are drawn on, and labelled with, the wrong arm.
#
# Nothing fitted is affected: the solid binned means come from the data, and
# every coefficient table indexes by explicit term name. This is the figure.
#
# Note that "alphabetically first" is not even fixed for a given cohort. For the
# melanoma arms it is "nonPD" under en_US collation and "PD" under C, so the
# same call can label the figure differently on two machines. The fixture below
# therefore decides at run time rather than hard-coding a level, and testthat
# runs under C.

condition_fixture <- function() {
  skip_if_not_installed("SpatialExperiment")
  path <- system.file("extdata", "mel_cosmx_subset.rds", package = "PACE")
  skip_if(path == "", "example dataset not installed")
  spe <- readRDS(path)
  spe <- spe[, spe$image %in% c("32158_19", "32176_5", "32177_4", "32178_7")]

  ## Make the model's reference the level that is NOT alphabetically first, so
  ## the model's order and the plot's old order are forced to disagree.
  arms <- sort(unique(as.character(spe$Responder)))
  expect_length(arms, 2L)
  spe$Responder <- stats::relevel(factor(spe$Responder), ref = arms[2])
  ## The premise of the whole test, asserted rather than assumed.
  expect_identical(levels(spe$Responder)[1], arms[2])
  expect_false(identical(levels(spe$Responder)[1], arms[1]))

  fit <- paceModel(spe, celltype_col = "cellType", condition_col = "Responder",
                   image_col = "image", n_iter = 2L, threads = 2L, verbose = FALSE)
  ## model.matrix() names the dummy for the NON-reference level, so the arm
  ## carrying the interaction is the alphabetically first one -- precisely the
  ## arm the old code handed the reference slope to.
  expect_identical(fit@params$resp_term, paste0("Responder", arms[1]))

  ## Pick the triple with the largest interaction, so the two arms cannot share
  ## a slope and the mapping is actually observable. A fixed choice is no good:
  ## the sparse-pair drop zeroes several pairs in a crop this small.
  u_matrix <- fit@fit$U
  marker <- paste0("::", fit@params$resp_term, ":")
  rows <- grep(marker, rownames(u_matrix), value = TRUE, fixed = TRUE)
  expect_gt(length(rows), 0L)
  block <- abs(u_matrix[rows, , drop = FALSE])
  at <- which(block == max(block), arr.ind = TRUE)[1, ]
  row <- rows[at[["row"]]]
  pair <- strsplit(sub(marker, "::", row, fixed = TRUE), "::", fixed = TRUE)[[1]]

  list(fit = fit, spe = spe, gene = colnames(u_matrix)[at[["col"]]],
       focal = pair[1], neighbour = pair[2], case_arm = arms[1], ref_arm = arms[2])
}

test_that("the dashed slopes follow the model's reference, not the alphabet", {
  skip_if_not_installed("ggplot2")
  fx <- condition_fixture()
  slopes <- PACE:::.pace_arm_slopes(fx$fit, fx$gene, fx$focal, fx$neighbour)
  ## Without a difference between the arms the assertions below would hold
  ## whichever way round the mapping went.
  expect_false(isTRUE(all.equal(slopes[["ref"]], slopes[["alt"]])))

  plot <- plotResponseCurve(fx$fit, fx$spe, fx$gene, fx$focal, fx$neighbour)

  ## The subtitle prints the case arm first. The case arm carries the
  ## interaction and so must take the alt slope; the reference takes ref.
  expect_identical(
    plot$labels$subtitle,
    sprintf("Solid = binned means +/- SE; dashed = PACE slope (%s %+.3f, %s %+.3f)",
            fx$case_arm, slopes[["alt"]], fx$ref_arm, slopes[["ref"]]))

  ## The labels are not the only thing that was wrong, so check the lines that
  ## were actually drawn. Found by linetype rather than by layer position.
  dashed <- vapply(plot$layers,
                   function(layer) identical(layer$aes_params$linetype, "dashed"),
                   logical(1))
  expect_identical(sum(dashed), 1L)
  drawn <- plot$layers[[which(dashed)]]$data
  drawn_slope <- function(arm) {
    rows <- drawn[as.character(drawn$arm) == arm, ]
    rows <- rows[order(rows$x), ]
    diff(rows$y) / diff(rows$x)
  }
  expect_equal(drawn_slope(fx$case_arm), slopes[["alt"]], tolerance = 1e-10)
  expect_equal(drawn_slope(fx$ref_arm), slopes[["ref"]], tolerance = 1e-10)
})

test_that("an arm that resp_term does not name is refused rather than guessed", {
  skip_if_not_installed("ggplot2")
  fx <- condition_fixture()

  broken <- fx$fit
  broken@params$resp_term <- "ResponderMISSING"
  expect_error(plotResponseCurve(broken, fx$spe, fx$gene, fx$focal, fx$neighbour),
               "not among the Responder levels present")

  broken@params$resp_term <- "SomethingElsePD"
  expect_error(plotResponseCurve(broken, fx$spe, fx$gene, fx$focal, fx$neighbour),
               "does not begin with condition_col")
})
