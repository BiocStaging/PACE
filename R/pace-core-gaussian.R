## ===========================================================================
## Gaussian (identity-link) PACE, for continuous intensities
## ===========================================================================
##
## DELIBERATELY SEPARATE from pace_fit_streaming(). The count path produced the
## numbers in a manuscript under review, so it carries no new branches for this:
## everything Gaussian lives here and the two are wired together only through
## the shared helpers they both call (pace_neighbour_kernel, pace_drop_sparse_k,
## pace_center_within_image, build_random_design_multi).
##
## Four things in the count orchestrator assume counts, and each is the reason a
## `family` flag there would not have been enough:
##
##   offset       pace_fit_streaming() sets offset = log(nCount), the library
##                size. On the identity link the offset is added to the mean
##                rather than multiplying it, and an intensity has no library
##                size, so it is ZERO here.
##   detection    the det_min re-filter keeps genes by detection RATE, which is
##                a count concept. Not applied.
##   contamination the Gaussian path is contamination-free BY CONSTRUCTION: the
##                technical contamination is removed from the intensities
##                upstream, against an ambient field, before they reach the fit.
##                There is no rho to estimate and no anchors to estimate it from.
##   mu           .pace_mu() rebuilds mu = exp(eta + log(nCount)). Here mu IS
##                eta + offset. This is the one that would have gone wrong
##                quietly: the fit would look healthy and the decomposition
##                built on mu would be wrong.
##
## NOTE: the random-effect formula construction below is DUPLICATED from
## pace_fit_streaming() (R/pace-core.R, "3. fixed effects + RE spec"). That is
## the price of keeping the count path untouched. If you change the formula
## logic in one, change it in the other.

#' Fit the PACE model to continuous intensities
#'
#' The Gaussian, identity-link form of [paceModel()], for continuous measurements
#' such as imaging mass cytometry protein intensities rather than counts.
#'
#' @details
#' The working response on the identity link is the intensity itself and the
#' working weight is the per-gene residual variance, so the mean converges in a
#' single inner solve and the outer iteration only moves the variance components.
#'
#' This path does not model technical contamination. Correct the intensities
#' first, against an ambient field built with [ambientField()], and pass the
#' corrected matrix here.
#'
#' @param object A [SpatialExperiment::SpatialExperiment].
#' @param celltype_col Column of `colData` holding the cell type labels.
#' @param image_col Column of `colData` identifying the image. Defaults to the
#'   single image when there is only one.
#' @param condition_col Optional column of `colData` giving a condition contrast.
#' @param assay_name The assay holding the intensities.
#' @param h_bio,h_tech Bandwidths of the biological and technical kernels.
#' @param kernel_per_image Build the kernels within each image separately.
#' @param image_re An optional second random-effect block over images.
#' @param types Cell types, in the order the fit should use them.
#' @param resp_term Optional response term passed to the model frame.
#' @param verbose Print the fitting trace.
#' @param ... Passed to `pace_fit_streaming_gaussian()`.
#'
#' @return A [PACEFit].
#' @rdname paceModelGaussian
#' @export
setMethod(
  "paceModelGaussian", "SpatialExperiment",
  function(object, celltype_col, image_col = NULL, condition_col = NULL,
           assay_name = "intensity", h_bio = 30, h_tech = 5,
           kernel_per_image = FALSE,
           image_re = c("none", "intercept", "slopes", "condition_slopes"),
           types = NULL, resp_term = NULL, verbose = TRUE, ...) {
    image_re <- match.arg(image_re)

    prep <- .pace_prepare(object, celltype_col, image_col, condition_col,
                          assay_name, resp_term)
    res <- pace_fit_streaming_gaussian(
      prep$Y, prep$df, types = types,
      celltype_col = celltype_col, image_col = prep$image_col,
      coord_cols = c("x", "y"), h_bio = h_bio, h_tech = h_tech,
      condition_col = condition_col, kernel_per_image = kernel_per_image,
      image_re = image_re, verbose = verbose, ...)

    methods::new(
      "PACEFit",
      fit       = res$fit,
      cellTypes = as.character(res$types),
      context   = list(df = res$df, X_fixed = res$X_fixed,
                       genes = colnames(res$fit$U)),
      params = list(h_bio = h_bio, h_tech = h_tech,
                    family = "gaussian",
                    ## Recorded so the readouts cannot mistake this for a count
                    ## fit: mu is eta + offset here, not exp(eta + offset).
                    contamination = "none", dispersion = "gaussian",
                    celltype_col = celltype_col, image_col = prep$image_col,
                    condition_col = condition_col, resp_term = prep$resp_term,
                    kernel_per_image = kernel_per_image, image_re = image_re,
                    assay_name = assay_name, edge_correct = FALSE))
  })


#' Gaussian streaming orchestrator
#'
#' Builds the neighbour kernels and the random-effect design for continuous
#' intensities, then fits them with the identity-link engine. The Gaussian
#' counterpart of `pace_fit_streaming()`; see the notes at the top of this file
#' for what differs and why they are kept apart.
#'
#' @param Y,df,types,celltype_col,image_col,coord_cols Intensities, the model
#'   frame, the cell types and the columns naming them.
#' @param h_bio,h_tech,eps Kernel bandwidths and the neighbour radius.
#' @param condition_col,kernel_per_image,image_re,drop_sparse_neff,within_image
#'   Design options, as `pace_fit_streaming()` takes them.
#' @param n_iter,threads,chunk_size,tau_shrinkage,early_stop_tol,min_iter,tau_max
#'   Fitting control, as `pace_fit_streaming()` takes them.
#' @param verbose Print the fitting trace.
#' @return A list with the fit, the model frame, the fixed-effect design, the
#'   intensities, the kernels and the cell types.
#' @keywords internal
pace_fit_streaming_gaussian <- function(Y, df, types = NULL,
                                        celltype_col, image_col,
                                        coord_cols = c("x", "y"),
                                        h_bio = 30, h_tech = 5, eps = NULL,
                                        condition_col = NULL,
                                        kernel_per_image = FALSE,
                                        image_re = c("none", "intercept",
                                                     "slopes", "condition_slopes"),
                                        drop_sparse_neff = 30,
                                        within_image = TRUE,
                                        n_iter = 32L, threads = 4L, chunk_size = 128L,
                                        tau_shrinkage = "adaptive",
                                        early_stop_tol = 2e-2, min_iter = 12L,
                                        tau_max = 100,
                                        verbose = TRUE) {
  image_re <- match.arg(image_re)
  if (is.null(eps)) eps <- 3 * h_bio
  .pace_check_positive(h_bio, "h_bio")
  .pace_check_positive(h_tech, "h_tech")
  .pace_check_positive(eps, "eps")
  has_cond <- !is.null(condition_col)
  if (image_re == "condition_slopes" && !has_cond)
    stop("`image_re = \"condition_slopes\"` needs a `condition_col`.", call. = FALSE)

  ## ---- 1. working frame ----
  ## No library size: see the note at the top of this file.
  df <- as.data.frame(df)
  celltype_raw <- as.character(df[[celltype_col]])
  if (is.null(types)) types <- sort(unique(celltype_raw))
  df$celltype <- factor(celltype_raw, levels = types)
  if (anyNA(df$celltype))
    stop(sum(is.na(df$celltype)), " cell(s) have a missing cell type or one not in `types`; ",
         "remove them or add their type to `types`.", call. = FALSE)
  df$imageID <- factor(as.character(df[[image_col]]))
  Y_sparse <- .pace_as_dgc(Y)
  coords <- as.matrix(df[, coord_cols])
  if (verbose)
    message(sprintf("pace_fit_streaming_gaussian: %d cells x %d markers; %d images",
                    nrow(Y_sparse), ncol(Y_sparse), nlevels(df$imageID)))

  ## ---- 2. neighbour kernels (shared with the count path) ----
  ker <- pace_neighbour_kernel(coords, df$celltype, types, h_bio, h_tech, eps,
                               image = df$imageID, per_image = kernel_per_image,
                               threads = threads)
  K_bio  <- ker$K_bio
  K_tech <- ker$K_tech
  if (drop_sparse_neff > 0)
    K_bio <- pace_drop_sparse_k(K_bio, df$celltype, types, drop_sparse_neff, verbose)
  if (within_image)
    K_bio <- pace_center_within_image(K_bio, df$celltype, df$imageID, types)
  for (tc in types) df[[tc]] <- K_bio[, tc]   ## re_specs formula reads these by name

  ## ---- 3. fixed effects + RE spec ----
  ## NOTE: duplicated from pace_fit_streaming(); keep the two in step.
  offset_vec <- numeric(nrow(df))             ## identity link: no offset
  X_fixed <- if (!has_cond) {
    matrix(1, nrow(df), 1, dimnames = list(NULL, "(Intercept)"))
  } else {
    stats::model.matrix(stats::as.formula(paste0("~ 1 + ", condition_col)), data = df)
  }

  types_rhs <- paste(types, collapse = " + ")
  celltype_formula <- if (has_cond)
    stats::as.formula(paste0("~ 1 + ", condition_col, " * (", types_rhs, ")"))
  else
    stats::as.formula(paste0("~ 1 + ", types_rhs))
  re_specs <- list(list(group_col = "celltype", formula = celltype_formula))

  if (image_re != "none") {
    ## A kernel column with no spread carries no slope to estimate, and
    ## standardising it would return NaN for every cell and drop every row from
    ## model.matrix(). Those columns are left out of the image-slope terms.
    standardised <- lapply(types, function(tc) pace_standardise_cpp(df[[tc]]))
    names(standardised) <- types
    varying <- vapply(standardised, function(z) is.finite(z$sd) && z$sd > 0, logical(1))
    if (image_re != "intercept" && any(!varying) && verbose)
      message("    image_re: no image slope for constant kernel column(s) ",
              paste(types[!varying], collapse = ", "))
    imgz <- sprintf("%s_imgz", types[varying])
    for (tc in types[varying]) df[[paste0(tc, "_imgz")]] <- standardised[[tc]]$values
    img_rhs <- switch(image_re,
      intercept        = "1",
      slopes           = paste(c("1", imgz), collapse = " + "),
      condition_slopes = paste(c("1", imgz, sprintf("%s:%s", condition_col, imgz)),
                               collapse = " + "))
    re_specs <- c(re_specs, list(list(
      group_col = "imageID",
      formula   = stats::as.formula(paste0("~ ", img_rhs)))))
  }

  re_design <- build_random_design_multi(df, re_specs)

  ## ---- 4. PQL fit (identity link) ----
  ## No ambient field, no anchors and no data-informed tau weights: the first two
  ## do not exist on this path, and the third is a count-variance heuristic.
  fit <- fit_pace_mvpql_streaming(
    Y = Y_sparse, X_fixed = X_fixed, df = df, re_specs = re_specs, re = re_design,
    offset_vec = offset_vec, data_informed_W = NULL,
    family = "gaussian",
    bleed_percell = FALSE,
    n_iter = as.integer(n_iter), tol = 5e-3,
    tau_shrinkage = tau_shrinkage,
    BPPARAM = BiocParallel::SerialParam(), n_threads = as.integer(threads),
    interior_precision = 1L, chunk_size = as.integer(chunk_size),
    early_stop_tol = early_stop_tol,
    min_iter = as.integer(min_iter),
    tau_max = tau_max,
    verbose = verbose)

  list(fit = fit, df = df, X_fixed = X_fixed, Y = Y_sparse,
       K_tech = K_tech, K_bio = K_bio, types = types)
}
