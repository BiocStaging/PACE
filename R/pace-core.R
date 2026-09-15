## pace_core.R -- cohort-agnostic PACE-MV method core (Stage 0 package extraction).
##
## A FAITHFUL lift of streaming/builders/build_bc_streaming.R: the computation and
## its order are identical, but every `Sys.getenv("R_*")` flag becomes an explicit
## function argument, and there is NO setwd / saveRDS / quit in the core. The
## solver itself (fit_pace_mvpql_streaming) and the reporting helpers
## (mvpql_to_results_multi, apply_mashr_shrinkage, mvpql_variance_decomposition_multi,
## build_random_design_multi, .compute_data_informed_weights) are assumed already
## sourced -- this file only repackages the BUILDER layer.
##
## Public surface (the future package API):
##   pace_neighbour_kernel()  -- Gaussian K_bio + exponential K_tech neighbour fields
##   pace_ambient_field()     -- sparse E^tech weight matrix W (streamed contamination)
##   pace_anchors()           -- homotypic-core negative-control anchors
##   pace_fit_streaming()     -- compose the above + streaming PQL fit  (THE method)
##   pace_shrink()            -- mash shrinkage of the neighbour slopes
##   pace_decompose()         -- per-gene variance decomposition (4-block)
##   pace_top_drivers()       -- per-pair MCSD driver tables
##
## Generic usage -- ANY spatial dataset, no per-cohort build script:
##   res <- pace_fit_streaming(Y, coldata,                 # counts + cell metadata
##                             celltype_col = "cellType",   # which columns hold the
##                             image_col    = "sample")     #   cell type + sample id
##   shr <- pace_shrink(res$fit, res$types)
##   dec <- pace_decompose(res$fit, res$df, res$Y, res$types, res$X_fixed)
##   ct <- as.character(res$df$celltype)
##   mu_means <- rowsum(res$fit$mu, ct) / as.vector(table(ct))   # rows named by cell type
##   drv <- pace_top_drivers(res$fit, shr, dec, res$types, mu_means)  # all focal x neighbour pairs
## Cell types, coordinates and image grouping are READ FROM THE DATA; the method
## defaults reproduce the manuscript recipe. (The Bioconductor wrapper paceFit(spe)
## reads these from a SpatialExperiment's conventions, collapsing the call to one
## argument.) Cohort-specific values (a fixed cell-type order, chosen pairs) are
## only ever needed to reproduce a *locked* fit -- see the reproduction test.

## ----------------------------------------------------------------------------
## Isotropic area-fraction edge correction.
## area_fraction(i) = area(disc(i, r) intersect image_rectangle) / (pi r^2),
## approximated by angular quadrature over `n_angles` directions.
## ----------------------------------------------------------------------------
pace_area_fraction <- function(coords, r,
                               xmin = NULL, xmax = NULL,
                               ymin = NULL, ymax = NULL,
                               n_angles = 1000L) {
  if (is.null(xmin)) xmin <- min(coords[, 1])
  if (is.null(xmax)) xmax <- max(coords[, 1])
  if (is.null(ymin)) ymin <- min(coords[, 2])
  if (is.null(ymax)) ymax <- max(coords[, 2])
  theta <- seq(0, 2 * pi, length.out = n_angles + 1L)[-1L]
  cos_t <- cos(theta)
  sin_t <- sin(theta)
  n <- nrow(coords)
  af <- numeric(n)
  norm_factor <- pi * r^2
  for (i in seq_len(n)) {
    x0 <- coords[i, 1]
    y0 <- coords[i, 2]
    d_right  <- ifelse(cos_t > 0, (xmax - x0) / cos_t, Inf)
    d_left   <- ifelse(cos_t < 0, (xmin - x0) / cos_t, Inf)
    d_top    <- ifelse(sin_t > 0, (ymax - y0) / sin_t, Inf)
    d_bottom <- ifelse(sin_t < 0, (ymin - y0) / sin_t, Inf)
    d_max <- pmin(d_right, d_left, d_top, d_bottom)
    r_eff <- pmin(r, d_max)
    af[i] <- (pi * sum(r_eff^2) / n_angles) / norm_factor
  }
  af
}

## ----------------------------------------------------------------------------
## Neighbour kernels: for every cell, the kernel-weighted abundance of each
## neighbour cell type within radius `eps`.
##   K_bio[i, c]  = sum_{j in type c, j != i} exp(-(d_ij / h_bio)^2)   (Gaussian)
##   K_tech[i, c] = sum_{j in type c, j != i} exp(-d_ij / h_tech)      (exponential)
## ----------------------------------------------------------------------------
pace_neighbour_kernel <- function(coords, celltype, types, h_bio, h_tech, eps,
                                  image = NULL, per_image = FALSE) {
  n <- nrow(coords)
  ct <- as.character(celltype)
  if (per_image) {
    ## images are separate samples (e.g. patients) -> no cross-image neighbours.
    stopifnot(!is.null(image))
    i_ok <- integer(0); jc_ok <- integer(0); d_ok <- numeric(0)
    for (im in unique(image)) {
      rows <- which(image == im)
      if (length(rows) < 2L) next
      fr <- dbscan::frNN(coords[rows, , drop = FALSE], eps = eps)
      i_loc <- rep.int(seq_along(rows), lengths(fr$id))
      j_loc <- unlist(fr$id,   use.names = FALSE)
      d_loc <- unlist(fr$dist, use.names = FALSE)
      jc <- match(ct[rows][j_loc], types)
      keep <- !is.na(jc)
      i_ok  <- c(i_ok,  rows[i_loc[keep]])
      jc_ok <- c(jc_ok, jc[keep])
      d_ok  <- c(d_ok,  d_loc[keep])
    }
  } else {
    ## one physical section -> global neighbours (cross-"slide" cells are real).
    fr <- dbscan::frNN(coords, eps = eps)
    i_vec <- rep.int(seq_len(n), lengths(fr$id))
    j_vec <- unlist(fr$id,   use.names = FALSE)
    d_vec <- unlist(fr$dist, use.names = FALSE)
    jc <- match(ct[j_vec], types)
    ok <- !is.na(jc)
    i_ok  <- i_vec[ok]
    d_ok  <- d_vec[ok]
    jc_ok <- jc[ok]
  }
  ## sparseMatrix SUMS duplicate (i, type) entries -> correct kernel accumulation
  K_tech <- as.matrix(Matrix::sparseMatrix(i = i_ok, j = jc_ok, x = exp(-d_ok / h_tech),
                                           dims = c(n, length(types))))
  K_bio <- as.matrix(Matrix::sparseMatrix(i = i_ok, j = jc_ok, x = exp(-d_ok^2 / h_bio^2),
                                          dims = c(n, length(types))))
  colnames(K_tech) <- paste0(types, "_near")
  colnames(K_bio)  <- types
  list(K_bio = K_bio, K_tech = K_tech)
}

## ----------------------------------------------------------------------------
## Drop (focal, neighbour) kernel columns with too little effective support.
## For each focal cell type, zero K_bio[focal cells, nb] when the centred
## column's effective sample size n_eff = sum(Kc^2) / max(Kc^2) < neff_min.
## ----------------------------------------------------------------------------
pace_drop_sparse_k <- function(K_bio, celltype, types, neff_min, verbose = TRUE) {
  n_dropped <- 0L
  for (focal in types) {
    cells_f <- which(celltype == focal)
    if (!length(cells_f)) next
    for (nb in types) {
      vals <- K_bio[cells_f, nb]
      centred <- vals - mean(vals, na.rm = TRUE)
      max_sq <- max(centred^2, na.rm = TRUE)
      n_eff <- if (max_sq > 0) sum(centred^2, na.rm = TRUE) / max_sq else 0
      if (n_eff < neff_min) {
        K_bio[cells_f, nb] <- 0
        n_dropped <- n_dropped + 1L
      }
    }
  }
  if (verbose)
    message(sprintf("    drop_sparse_k: zeroed %d (focal, neighbour) pairs with n_eff < %g",
                    n_dropped, neff_min))
  K_bio
}

## ----------------------------------------------------------------------------
## Within-(image, celltype) centring of each neighbour-kernel column, so the
## random slopes are estimated on the within-group deviation only.
## ----------------------------------------------------------------------------
pace_center_within_image <- function(K_bio, celltype, image, types) {
  for (tc in types) {
    col <- K_bio[, tc]
    K_bio[, tc] <- col - stats::ave(col, image, celltype, FUN = mean)
  }
  K_bio
}

## ----------------------------------------------------------------------------
## Homotypic-core negative-control anchors (de-circularised contamination refs).
## A cell type's clean identity profile is estimated from its spatially-isolated
## (>= homo_frac same-type-neighbour) cells; a gene anchors type X when another
## type clearly owns it (owner_mean > owner_thresh) and X's core level is < a
## small fraction (core_thresh) of the owner's.
## Returns the memory-light form: an (n_types x G) 0/1 mask + a length-n celltype
## index (1..n_types), expanded per cell inside the solver.
## ----------------------------------------------------------------------------
pace_anchors <- function(coords, Y, celltype, image, types,
                         homo_frac = 0.5, owner_thresh = 0.1, core_thresh = 0.1,
                         verbose = TRUE) {
  n <- nrow(Y)
  type_means <- t(vapply(types,
                         function(tt) colMeans(Y[celltype == tt, , drop = FALSE]),
                         numeric(ncol(Y))))
  rownames(type_means) <- types
  colnames(type_means) <- colnames(Y)

  ## same-type-neighbour fraction within 30 um, per image
  same_frac <- rep(NA_real_, n)
  for (s in unique(image)) {
    si <- which(image == s)
    if (length(si) < 50) next
    ctl <- celltype[si]
    fr  <- dbscan::frNN(coords[si, , drop = FALSE], eps = 30)
    for (k in seq_along(si)) {
      idk <- fr$id[[k]]
      same_frac[si[k]] <- if (length(idk)) mean(ctl[idk] == ctl[k]) else 1
    }
  }
  core <- which(same_frac >= homo_frac)

  ## clean profile from core cells (raw type-mean fallback when a type is sparse)
  core_means <- type_means
  for (X in types) {
    idx <- core[celltype[core] == X]
    if (length(idx) >= 20) core_means[X, ] <- colMeans(Y[idx, , drop = FALSE])
  }
  owner_mean <- apply(core_means, 2, max)
  owner_t    <- types[apply(core_means, 2, which.max)]

  mask <- matrix(0, length(types), ncol(Y), dimnames = list(types, colnames(Y)))
  n_anchor <- integer(0)
  for (ti in seq_along(types)) {
    X <- types[ti]
    is_anchor <- owner_t != X &
                 owner_mean > owner_thresh &
                 (core_means[X, ] / pmax(owner_mean, 1e-9) < core_thresh)
    mask[ti, ] <- as.numeric(is_anchor)
    n_anchor <- c(n_anchor, sum(is_anchor))
  }
  if (verbose)
    message(sprintf("    anchors: homotypic-core cells %d/%d; anchors/type median=%d range=[%d,%d]",
                    length(core), n, as.integer(stats::median(n_anchor)),
                    min(n_anchor), max(n_anchor)))
  list(mask = mask, idx = as.integer(celltype))
}

## ----------------------------------------------------------------------------
## Sparse E^tech ambient weight matrix W (n x n).
##   W[i, j] = exp(-d_ij / h_tech) / area_fraction(i)   for cross-celltype
##             neighbours j within rad = 3 * h_tech of i (self excluded).
## Streaming identity: W %*% Y[, g] == dense E_tech[, g]. Validated on a few
## random genes before the fit when `validate = TRUE`.
## ----------------------------------------------------------------------------
pace_ambient_field <- function(coords, Y, celltype, image, types, h_tech,
                               edge_correct = TRUE, validate = TRUE, verbose = TRUE) {
  n <- nrow(Y)
  rad <- 3 * h_tech
  by_image <- split(seq_len(n), image)
  ct_all <- as.character(celltype)

  ## per-image edge area-fraction at r = rad
  af <- rep(1, n)
  if (edge_correct) {
    for (si in seq_along(by_image)) {
      rows <- by_image[[si]]
      if (!length(rows)) next
      ci <- coords[rows, , drop = FALSE]
      af[rows] <- pace_area_fraction(ci, rad,
                                     min(ci[, 1]), max(ci[, 1]),
                                     min(ci[, 2]), max(ci[, 2]))
    }
  }

  ## (i, j, w) triplets for the cross-celltype neighbour weights
  ii_all <- vector("list", length(by_image))
  jj_all <- vector("list", length(by_image))
  ww_all <- vector("list", length(by_image))
  for (si in seq_along(by_image)) {
    rows <- by_image[[si]]
    if (length(rows) < 2) next
    ci <- coords[rows, , drop = FALSE]
    nn <- dbscan::frNN(ci, eps = rad)
    ct_im <- ct_all[rows]
    i_loc <- rep.int(seq_along(rows), lengths(nn$id))
    j_loc <- unlist(nn$id,   use.names = FALSE)
    d_loc <- unlist(nn$dist, use.names = FALSE)
    keep  <- ct_im[j_loc] != ct_im[i_loc]
    if (!any(keep)) next
    i_loc <- i_loc[keep]
    j_loc <- j_loc[keep]
    d_loc <- d_loc[keep]
    i_glb <- rows[i_loc]
    j_glb <- rows[j_loc]
    w_glb <- exp(-d_loc / h_tech) / af[i_glb]   ## edge correction folds into row i
    ii_all[[si]] <- i_glb
    jj_all[[si]] <- j_glb
    ww_all[[si]] <- w_glb
  }
  ii <- unlist(ii_all, use.names = FALSE)
  jj <- unlist(jj_all, use.names = FALSE)
  ww <- unlist(ww_all, use.names = FALSE)
  W <- Matrix::sparseMatrix(i = ii, j = jj, x = ww, dims = c(n, n))
  W <- methods::as(W, "CsparseMatrix")

  ## cheap correctness gate: W %*% Y == dense E_tech on a few spot-check genes.
  ## Genes are chosen deterministically (evenly spaced) so the check does not
  ## touch the global RNG state.
  if (validate) {
    spot_g <- unique(round(seq(1, ncol(Y), length.out = min(5L, ncol(Y)))))
    E_spot <- matrix(0, n, length(spot_g))
    for (si in seq_along(by_image)) {
      rows <- by_image[[si]]
      if (length(rows) < 2) next
      ci <- coords[rows, , drop = FALSE]
      nn <- dbscan::frNN(ci, eps = rad)
      Y_im  <- Y[rows, spot_g, drop = FALSE]
      ct_im <- ct_all[rows]
      for (i in seq_along(rows)) {
        nbr <- nn$id[[i]]
        if (!length(nbr)) next
        dl <- nn$dist[[i]]
        kk <- ct_im[nbr] != ct_im[i]
        if (!any(kk)) next
        w <- exp(-dl[kk] / h_tech)
        E_spot[rows[i], ] <- as.numeric(crossprod(w, Y_im[nbr[kk], , drop = FALSE]))
      }
      E_spot[rows, ] <- E_spot[rows, , drop = FALSE] / af[rows]
    }
    a_spot <- as.matrix(W %*% Y[, spot_g, drop = FALSE])
    max_id <- max(abs(a_spot - E_spot))
    if (verbose)
      message(sprintf("    [W identity check] max|W%%*%%Y - dense E_tech| over %d genes = %.3e",
                      length(spot_g), max_id))
    if (max_id > 1e-8)
      stop("pace_ambient_field: W identity check FAILED (max abs diff ", max_id, ")")
  }

  list(W = W, image_idx = as.integer(image), n_images = nlevels(image))
}

## ----------------------------------------------------------------------------
## Fitted means and the contamination log-offset, rebuilt rather than stored.
##
## mu is a deterministic function of what the fit already holds: on the log
## scale eta = X B + Z U has rank at most p + q, so B, U, X_fixed and re_meta$Z
## are the factored form of it, and the dense n x G matrix is the expanded copy.
## Only the ambient field has to be recomputed, from the same counts and the
## same settings the fit used. Reconstruction is exact, not approximate.
##
## technical_offset_mat = log1p(mu_spill / mu_bio) must be rebuilt alongside mu:
## it is what gates the spillover block of the decomposition (`use_bleed`), so a
## fit missing it does not fail, it silently reports no spillover.
##
## The settings must come from the fit, never from the defaults of the exported
## ambientField(): a fit with another h_tech, image grouping, cell-type order or
## edge correction would otherwise be handed a different field and would return
## plausible, wrong numbers.
## ----------------------------------------------------------------------------
.pace_mu_parts <- function(object, spe) {
  f <- object@fit
  ## The two travel together: both solvers return them from the same pass, so
  ## either both are present or neither is.
  if (!is.null(f$mu) && !is.null(f$technical_offset_mat))
    return(list(mu = f$mu, technical_offset_mat = f$technical_offset_mat))

  inputs <- .pace_mu_inputs(object, spe)
  inputs$Y <- as.matrix(inputs$Y)
  block <- .pace_mu_block(object, seq_len(nrow(inputs$Y)), inputs)

  ## contamination = "none" carries no ambient term: mu = mu_bio, no spillover.
  if (is.null(block$mu_spill))
    return(list(mu = block$mu_bio,
                technical_offset_mat = matrix(0, nrow(block$mu_bio), ncol(block$mu_bio),
                                              dimnames = dimnames(block$mu_bio))))

  list(mu                   = pmax(block$mu_bio + block$mu_spill, 1e-6),
       technical_offset_mat = log1p(block$mu_spill / block$mu_bio))
}

## Mean fitted mean of every gene within each cell type (cell types x genes),
## computed chunk by chunk so the n x G `mu` is never held whole. The driver
## scores need only these means, and at a million cells each n x G matrix runs
## to gigabytes. A fit that stores `mu` is averaged directly.
.pace_mu_means_by_celltype <- function(object, spe = NULL, chunk_size = 10000L) {
  f <- object@fit
  types <- object@cellTypes

  ## the cell-type membership pace_top_drivers() uses: the intercept column
  Z_intercepts <- f$re_meta$Z[, paste0(types, "::(Intercept)"), drop = FALSE]
  cell_type <- max.col(as.matrix(Z_intercepts != 0), ties.method = "first")
  cell_type[Matrix::rowSums(Z_intercepts != 0) == 0] <- NA_integer_
  cell_type <- factor(types[cell_type], levels = types)
  n_cells <- nrow(Z_intercepts)

  sums <- NULL
  add_block <- function(mu_block, rows) {
    typed <- !is.na(cell_type[rows])
    block_sums <- rowsum(mu_block[typed, , drop = FALSE], cell_type[rows][typed],
                         reorder = TRUE)
    full <- matrix(0, length(types), ncol(mu_block),
                   dimnames = list(types, colnames(mu_block)))
    full[rownames(block_sums), ] <- block_sums
    if (is.null(sums)) full else sums + full
  }

  if (!is.null(f$mu)) {
    sums <- add_block(as.matrix(f$mu), seq_len(n_cells))
  } else {
    if (is.null(spe))
      stop("this fit was saved without its fitted means; pass the ",
           "SpatialExperiment it was fitted on.", call. = FALSE)
    inputs <- .pace_mu_inputs(object, spe)
    for (start in seq(1L, n_cells, by = chunk_size)) {
      rows <- start:min(start + chunk_size - 1L, n_cells)
      block <- .pace_mu_block(object, rows, inputs)
      mu_block <- if (is.null(block$mu_spill)) block$mu_bio
                  else pmax(block$mu_bio + block$mu_spill, 1e-6)
      sums <- add_block(mu_block, rows)
    }
  }
  sums / as.vector(table(cell_type))
}

## What rebuilding mu needs beyond the fit: the counts (cells x fitted genes,
## kept sparse if the assay is) and the ambient field, or NULL for a fit
## without contamination.
.pace_mu_inputs <- function(object, spe) {
  f  <- object@fit
  df <- object@context$df
  p  <- object@params
  counts <- SummarizedExperiment::assay(spe, p$assay_name)
  if (ncol(counts) != nrow(df))
    stop("`spe` has ", ncol(counts), " cells but the fit has ", nrow(df),
         "; pass the same object used for paceModel().", call. = FALSE)
  Y <- Matrix::t(counts[object@context$genes, , drop = FALSE])

  if (!identical(p$contamination, "percell_hc") || is.null(f$percell_bleed_rho))
    return(list(Y = Y, amb = NULL))

  if (is.null(p$edge_correct))
    stop("this fit predates the recording of `edge_correct` and its ambient ",
         "field cannot be rebuilt faithfully; refit, or decompose a fit that ",
         "retains `mu`.", call. = FALSE)

  amb <- pace_ambient_field(SpatialExperiment::spatialCoords(spe), Y,
                            df$celltype, df$imageID, object@cellTypes,
                            h_tech = p$h_tech, edge_correct = p$edge_correct,
                            validate = FALSE, verbose = FALSE)
  list(Y = Y, amb = amb)
}

## The two parts of the fitted mean for a block of cells (rows x genes):
## mu_bio, and mu_spill (NULL without contamination). The one place the
## rebuild formula lives, shared by the whole-matrix and chunked callers.
.pace_mu_block <- function(object, rows, inputs) {
  f  <- object@fit
  df <- object@context$df

  ## mu_bio = exp(eta + offset), offset = log library size (pace_fit_streaming).
  eta    <- as.matrix(object@context$X_fixed[rows, , drop = FALSE] %*% f$B) +
            as.matrix(f$re_meta$Z[rows, , drop = FALSE] %*% f$U)
  mu_bio <- pmax(exp(eta + log(df$nCount[rows])), 1e-6)
  if (is.null(inputs$amb))
    return(list(mu_bio = mu_bio, mu_spill = NULL))

  ## Same expressions and floors as the solver's final pass.
  ambient  <- as.matrix(inputs$amb$W[rows, , drop = FALSE] %*% inputs$Y)
  mu_spill <- pmax(ambient * f$percell_bleed_rho[rows], 0)
  list(mu_bio = mu_bio, mu_spill = mu_spill)
}

## ----------------------------------------------------------------------------
## THE method: compose kernels + ambient + anchors + data-informed tau and run
## the streaming PQL fit. All knobs that were `R_*` env flags in the builder are
## explicit arguments here. Returns the fit plus the working frame / X_fixed / Y
## needed by the reporting functions below.
## ----------------------------------------------------------------------------
pace_fit_streaming <- function(Y, df, types = NULL,
                               celltype_col, image_col, coord_cols = c("x", "y"),
                               h_bio = 30, h_tech = 5, eps = NULL,
                               contamination = c("percell_hc", "none"),
                               dispersion = c("nb1", "nb2"),              ## NB1 is canonical; matches paceModel()
                               condition_col = NULL,                      ## disease/condition column (optional)
                               kernel_per_image = FALSE,                  ## TRUE when images = separate samples
                               image_re = c("none", "intercept",          ## second RE block over images
                                            "slopes", "condition_slopes"),
                               drop_sparse_neff = 30,
                               within_image = TRUE,
                               edge_correct = TRUE,
                               data_informed_tau = TRUE,
                               det_min = 0.05, homo_frac = 0.5,
                               n_iter = 32L, threads = 4L, chunk_size = 128L,
                               tau_shrinkage = "adaptive",
                               alpha_warmup = 6, early_stop_tol = 2e-2, min_iter = 12L,
                               fuse = FALSE, return_mu = TRUE,
                               ## Upper bound on the variance components, passed to the fitter.
                               ## A binding cap means the term is identified only by the ridge.
                               tau_max = 100,
                               verbose = TRUE) {
  contamination <- match.arg(contamination)
  dispersion    <- match.arg(dispersion)
  image_re      <- match.arg(image_re)
  if (is.null(eps)) eps <- 3 * h_bio
  use_etech <- contamination == "percell_hc"   ## E^tech ambient drives spillover
  has_cond  <- !is.null(condition_col)

  ## ---- 1. working frame: raw labels -> factors, library size ----
  df <- as.data.frame(df)
  celltype_raw <- as.character(df[[celltype_col]])
  ## default: every observed cell type (a fixed order is only needed to reproduce
  ## a locked fit, in which case the caller passes `types` explicitly).
  if (is.null(types)) types <- sort(unique(celltype_raw))
  df$celltype <- factor(celltype_raw, levels = types)
  df$imageID  <- factor(as.character(df[[image_col]]))
  ## optional max-per-celltype detection re-filter (no-op at det_min = 0.05)
  if (det_min > 0.05 + 1e-6) {
    ct_chr <- as.character(df$celltype)
    ct_idx <- lapply(types, function(c) which(ct_chr == c))
    max_det <- vapply(seq_len(ncol(Y)), function(gi)
      max(vapply(ct_idx, function(idx)
        if (length(idx) == 0L) 0 else mean(Y[idx, gi] > 0), numeric(1))), numeric(1))
    Y <- Y[, which(max_det >= det_min), drop = FALSE]
  }
  Y <- as.matrix(Y)
  df$nCount <- rowSums(Y)
  coords <- as.matrix(df[, coord_cols])
  if (verbose)
    message(sprintf("pace_fit_streaming: %d cells x %d genes; %d images",
                    nrow(Y), ncol(Y), nlevels(df$imageID)))

  ## ---- 2. neighbour kernels (+ sparse-pair drop + within-image centring) ----
  ker <- pace_neighbour_kernel(coords, df$celltype, types, h_bio, h_tech, eps,
                               image = df$imageID, per_image = kernel_per_image)
  K_bio  <- ker$K_bio
  K_tech <- ker$K_tech
  if (drop_sparse_neff > 0)
    K_bio <- pace_drop_sparse_k(K_bio, df$celltype, types, drop_sparse_neff, verbose)
  if (within_image)
    K_bio <- pace_center_within_image(K_bio, df$celltype, df$imageID, types)
  for (tc in types) df[[tc]] <- K_bio[, tc]   ## re_specs formula reads these by name

  ## ---- 3. fixed effects + RE spec ----
  ## E^tech path: X_fixed is intercept-only, or + a condition main effect when a
  ## disease contrast is supplied. The celltype RE block carries the neighbour
  ## slopes, crossed with the condition when present.
  offset_vec <- log(df$nCount)
  X_fixed <- if (use_etech && !has_cond) {
    matrix(1, nrow(df), 1, dimnames = list(NULL, "(Intercept)"))
  } else if (use_etech && has_cond) {
    stats::model.matrix(stats::as.formula(paste0("~ 1 + ", condition_col)), data = df)
  } else {
    cbind(`(Intercept)` = 1, K_tech)
  }
  types_rhs <- paste(types, collapse = " + ")
  celltype_formula <- if (has_cond)
    stats::as.formula(paste0("~ 1 + ", condition_col, " * (", types_rhs, ")"))
  else
    stats::as.formula(paste0("~ 1 + ", types_rhs))
  re_specs <- list(list(group_col = "celltype", formula = celltype_formula))

  ## optional SECOND RE block over images (e.g. patient-level neighbour slopes).
  ## Uses standardised kernel columns <type>_imgz; condition_slopes adds the
  ## condition interaction (the PAT_DZ design).
  if (image_re != "none") {
    imgz <- paste0(types, "_imgz")
    for (tc in types) df[[paste0(tc, "_imgz")]] <- as.numeric(scale(df[[tc]]))
    img_rhs <- switch(image_re,
      intercept        = "1",
      slopes           = paste(c("1", imgz), collapse = " + "),
      condition_slopes = paste(c("1", imgz, paste0(condition_col, ":", imgz)), collapse = " + "))
    re_specs <- c(re_specs, list(list(
      group_col = "imageID",
      formula   = stats::as.formula(paste0("~ ", img_rhs)))))
  }

  ## ---- 4. anchors + sparse ambient field ----
  anchors <- if (contamination == "percell_hc")
    pace_anchors(coords, Y, df$celltype, df$imageID, types, homo_frac, verbose = verbose)
  else NULL
  ambient_W <- NULL
  ambient_image_idx <- NULL
  ambient_n_images <- 0L
  if (use_etech) {
    amb <- pace_ambient_field(coords, Y, df$celltype, df$imageID, types, h_tech,
                              edge_correct = edge_correct, verbose = verbose)
    ambient_W <- amb$W
    ambient_image_idx <- amb$image_idx
    ambient_n_images <- amb$n_images
  }

  ## ---- 5. data-informed tau weights ----
  data_informed_W <- NULL
  if (data_informed_tau) {
    re_tmp <- build_random_design_multi(df, re_specs)
    data_informed_W <- .compute_data_informed_weights(
      re = re_tmp, Y = Y, df = df,
      focals = re_tmp$blocks[[1]]$group_levels, TYPES = types,
      celltype_col = "celltype", verbose = verbose)
    rm(re_tmp)
  }

  ## ---- 6. PQL fit ----
  ## The streaming solver only implements the per-cell contamination (additive)
  ## path: it streams the ambient field as W %*% Y. With contamination = "none"
  ## there is no ambient field to stream, so we fall back to the dense oracle
  ## solver (the same solver the streaming path is byte-identical against;
  ## feasible for targeted panels). The percell_hc path is unchanged.
  if (contamination == "percell_hc") {
    Y_sparse <- methods::as(methods::as(methods::as(Y, "dMatrix"), "generalMatrix"),
                            "CsparseMatrix")
    if (!is.null(colnames(Y))) colnames(Y_sparse) <- colnames(Y)
    fit <- fit_pace_mvpql_streaming(
      Y = Y_sparse, X_fixed = X_fixed, df = df, re_specs = re_specs,
      offset_vec = offset_vec, data_informed_W = data_informed_W,
      ambient_W = ambient_W, ambient_image_idx = ambient_image_idx,
      ambient_n_images = ambient_n_images,
      bleed_percell = TRUE,
      percell_anchor_mask = anchors$mask, percell_anchor_idx = anchors$idx,
      n_iter = as.integer(n_iter), tol = 5e-3,
      disp_model = dispersion,
      tau_shrinkage = tau_shrinkage,
      BPPARAM = BiocParallel::SerialParam(), n_threads = as.integer(threads),
      interior_precision = 1L, chunk_size = as.integer(chunk_size),
      alpha_warmup = alpha_warmup, early_stop_tol = early_stop_tol,
      min_iter = as.integer(min_iter), fuse_rho = fuse,
      tau_max = tau_max,
      return_mu = return_mu, verbose = verbose)
  } else {
    ## contamination == "none": no ambient field; dense oracle (bleed_percell = FALSE).
    fit <- fit_pace_mvpql_joint_multi(
      Y = Y, X_fixed = X_fixed, df = df, re_specs = re_specs,
      offset_vec = offset_vec, data_informed_W = data_informed_W,
      n_iter = as.integer(n_iter), tol = 5e-3,
      disp_model = dispersion, tau_shrinkage = tau_shrinkage,
      BPPARAM = BiocParallel::SerialParam(), n_threads = as.integer(threads),
      interior_precision = 1L, chunk_size = as.integer(chunk_size),
      tau_max = tau_max, verbose = verbose)
  }

  list(fit = fit, df = df, X_fixed = X_fixed, Y = Y,
       K_tech = K_tech, K_bio = K_bio, types = types,
       ## Returned so the fit can record it: it is the one ambient-field input
       ## not otherwise recoverable from the fit, and .pace_mu() needs it to
       ## rebuild the same field.
       edge_correct = edge_correct)
}

## ----------------------------------------------------------------------------
## Reporting: mash shrinkage of the neighbour slopes.
## ----------------------------------------------------------------------------
pace_shrink <- function(fit, types, resp_term = NULL, ...) {
  results_mv <- mvpql_to_results_multi(fit, keep_block = "celltype")
  apply_mashr_shrinkage(results = results_mv, focals = types,
                        neighbours = types, resp_term = resp_term, ...)
}

## ----------------------------------------------------------------------------
## Reporting: per-gene variance decomposition with the 4-block percentage view.
## ----------------------------------------------------------------------------
pace_decompose <- function(fit, df, Y, types, X_fixed, resp_term = NULL) {
  dec <- mvpql_variance_decomposition_multi(
    fit = fit, df = df, Y = Y, vars = types,
    X_fixed = X_fixed, resp_term = resp_term, focal_levels = types)
  g5 <- dec$gene_focal_5block
  total4 <- with(g5, celltype_offset_sq + V_state_baseline + V_state_responder +
                     V_spill + V_disp)
  dec$gene_focal_4block <- tibble::tibble(
    focal = g5$focal,
    gene  = g5$gene,
    `Cell type %`     = 100 * g5$celltype_offset_sq / pmax(total4, 1e-12),
    `Spatial state %` = 100 * (g5$V_state_baseline + g5$V_state_responder) / pmax(total4, 1e-12),
    `Spillover %`     = 100 * g5$V_spill / pmax(total4, 1e-12),
    `Residual %`      = 100 * g5$V_disp  / pmax(total4, 1e-12),
    celltype_offset_sq = g5$celltype_offset_sq,
    V_state = g5$V_state_baseline + g5$V_state_responder,
    V_spill = g5$V_spill,
    V_disp  = g5$V_disp,
    Total   = total4,
    spec       = g5$spec,
    focal_mean = g5$focal_mean)
  dec
}

## ----------------------------------------------------------------------------
## Reporting: per-pair MCSD driver tables.
##   MCSD = b_shrunk^2 * spec^2 * focal_mean, filtered lfsr < 0.05, ranked desc.
## `pairs` is a list of c(focal, neighbour); resp_term = NULL gives term == nb.
## ----------------------------------------------------------------------------
## `mu_means` is the cell types x genes matrix of mean fitted means
## (.pace_mu_means_by_celltype()), so the n x G `mu` need not be held.
pace_top_drivers <- function(fit, shrunken_long, dec, types, mu_means, pairs = NULL,
                             resp_term = NULL, resp_dummy = NULL) {
  ## default: every ordered focal != neighbour pair (a chosen subset is only a
  ## reporting convenience, not part of the method).
  if (is.null(pairs)) {
    pairs <- list()
    for (fc in types) for (nc in types) if (fc != nc) pairs <- c(pairs, list(c(fc, nc)))
  }
  g5 <- dec$gene_focal_5block
  Z_re <- fit$re_meta$Z
  cells_by_ct <- lapply(types, function(c) which(Z_re[, paste0(c, "::(Intercept)")] != 0))
  names(cells_by_ct) <- types
  alpha_g <- pmax(fit$alpha, 0)
  gene_names_fit <- colnames(fit$U)

  out <- list()
  for (p in pairs) {
    fc <- p[1]
    nc <- p[2]
    pk <- paste(fc, nc, sep = "_")
    target_term <- if (is.null(resp_term)) nc else paste0(resp_term, ":", nc)
    s <- shrunken_long |>
      dplyr::filter(focal == fc, neighbour == nc, term == target_term) |>
      dplyr::distinct(gene, .keep_all = TRUE)
    fm <- g5 |>
      dplyr::filter(focal == fc) |>
      dplyr::select(gene, spec, focal_mean) |>
      dplyr::distinct(gene, .keep_all = TRUE)
    if (!nrow(s) || !nrow(fm)) next

    cells_c <- cells_by_ct[[fc]]
    col_N <- paste0(fc, "::", nc)
    if (!col_N %in% colnames(Z_re)) {
      out[[pk]] <- list(scores = s[0, ], status = "dropped (n_eff)",
                        expected_false_sign = 0, false_sign_rate = NA_real_)
      next
    }
    N_t <- as.numeric(Z_re[cells_c, col_N])
    var_N <- stats::var(N_t, na.rm = TRUE)
    mu_bar_per_gene <- mu_means[fc, ]
    names(alpha_g) <- gene_names_fit
    names(mu_bar_per_gene) <- gene_names_fit

    base <- s |>
      dplyr::inner_join(fm, by = "gene") |>
      dplyr::mutate(
        MCSD    = (estimate_shrunk^2) * (spec^2) * pmax(focal_mean, 0),
        MCSD4   = (estimate_shrunk^2) * (spec^4) * pmax(focal_mean, 0),
        mu_bar  = mu_bar_per_gene[gene],
        alpha   = alpha_g[gene],
        V_resid = log(1 + (1 + pmax(alpha, 0)) / pmax(mu_bar, 1e-6)))
    if (is.null(resp_term)) {
      ## baseline neighbour effect only (no condition).
      res_all <- base |>
        dplyr::mutate(V_S = estimate_shrunk^2 * var_N,
                      V_total = V_S + V_resid,
                      R2_S = V_S / pmax(V_total, 1e-12))
      res <- res_all |>
        dplyr::filter(lfsr < 0.05) |>
        dplyr::arrange(dplyr::desc(MCSD)) |>
        dplyr::mutate(rank = dplyr::row_number()) |>
        dplyr::rename(b_clean = estimate_shrunk) |>
        dplyr::select(rank, gene, MCSD, MCSD4, b_clean, spec, focal_mean,
                      R2_S, mu_bar, alpha, V_S, V_resid, V_total, lfsr, sd_shrunk)
    } else {
      ## condition x spatial: baseline slope V_S (from the raw BLUP u) plus the
      ## condition-interaction slope V_RxS (from the shrunken estimate).
      col_R <- match(paste0(fc, "::", resp_term), colnames(Z_re))
      R_c <- if (is.na(col_R)) {
        if (is.null(resp_dummy)) rep(0, length(cells_c)) else resp_dummy[cells_c]
      } else as.numeric(Z_re[cells_c, col_R])
      var_RN <- stats::var(R_c * N_t, na.rm = TRUE)
      row_u <- match(col_N, rownames(fit$U))
      u_vec <- if (is.na(row_u))
        setNames(rep(0, length(gene_names_fit)), gene_names_fit)
      else setNames(as.numeric(fit$U[row_u, ]), gene_names_fit)
      res_all <- base |>
        dplyr::mutate(u_raw   = u_vec[gene],
                      V_S     = u_raw^2 * var_N,
                      V_RxS   = estimate_shrunk^2 * var_RN,
                      V_total = V_S + V_RxS + V_resid,
                      R2_S    = V_S   / pmax(V_total, 1e-12),
                      R2_RxS  = V_RxS / pmax(V_total, 1e-12))
      res <- res_all |>
        dplyr::filter(lfsr < 0.05) |>
        dplyr::arrange(dplyr::desc(MCSD)) |>
        dplyr::mutate(rank = dplyr::row_number()) |>
        dplyr::rename(b_clean = estimate_shrunk) |>
        dplyr::select(rank, gene, MCSD, MCSD4, b_clean, u_raw, spec, focal_mean,
                      R2_S, R2_RxS, mu_bar, alpha, V_S, V_RxS, V_resid, V_total,
                      lfsr, sd_shrunk)
    }
    status <- if (nrow(res) >= 3) "significant" else "honestly null"
    ## how many of this pair's calls are expected to have the wrong sign
    fsr <- expected_false_sign(res$lfsr)
    out[[pk]] <- list(scores = res, status = status,
                      expected_false_sign = fsr$expected_false_sign,
                      false_sign_rate = fsr$false_sign_rate)
  }
  out
}
