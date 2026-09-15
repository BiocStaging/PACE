# Reference implementations of the neighbourhood, edge-correction, ambient-field
# and anchor computations, copied verbatim from R/pace-core.R at f55d976 (the
# pure-R versions the compiled core replaces). Tests compare the compiled code
# against these, so the comparisons need no external oracle.

reference_area_fraction <- function(coords, r,
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

reference_neighbour_kernel <- function(coords, celltype, types, h_bio, h_tech, eps,
                                       image = NULL, per_image = FALSE) {
  n <- nrow(coords)
  ct <- as.character(celltype)
  if (per_image) {
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
  K_tech <- as.matrix(Matrix::sparseMatrix(i = i_ok, j = jc_ok, x = exp(-d_ok / h_tech),
                                           dims = c(n, length(types))))
  K_bio <- as.matrix(Matrix::sparseMatrix(i = i_ok, j = jc_ok, x = exp(-d_ok^2 / h_bio^2),
                                          dims = c(n, length(types))))
  colnames(K_tech) <- paste0(types, "_near")
  colnames(K_bio)  <- types
  list(K_bio = K_bio, K_tech = K_tech)
}

reference_anchors <- function(coords, Y, celltype, image, types,
                              homo_frac = 0.5, owner_thresh = 0.1, core_thresh = 0.1) {
  n <- nrow(Y)
  type_means <- t(vapply(types,
                         function(tt) colMeans(Y[celltype == tt, , drop = FALSE]),
                         numeric(ncol(Y))))
  rownames(type_means) <- types
  colnames(type_means) <- colnames(Y)
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
  core_means <- type_means
  for (X in types) {
    idx <- core[celltype[core] == X]
    if (length(idx) >= 20) core_means[X, ] <- colMeans(Y[idx, , drop = FALSE])
  }
  owner_mean <- apply(core_means, 2, max)
  owner_t    <- types[apply(core_means, 2, which.max)]
  mask <- matrix(0, length(types), ncol(Y), dimnames = list(types, colnames(Y)))
  for (ti in seq_along(types)) {
    X <- types[ti]
    is_anchor <- owner_t != X &
                 owner_mean > owner_thresh &
                 (core_means[X, ] / pmax(owner_mean, 1e-9) < core_thresh)
    mask[ti, ] <- as.numeric(is_anchor)
  }
  list(mask = mask, idx = as.integer(celltype), same_frac = same_frac)
}

reference_ambient_field <- function(coords, Y, celltype, image, types, h_tech,
                                    edge_correct = TRUE) {
  n <- nrow(Y)
  rad <- 3 * h_tech
  by_image <- split(seq_len(n), image)
  ct_all <- as.character(celltype)
  af <- rep(1, n)
  if (edge_correct) {
    for (si in seq_along(by_image)) {
      rows <- by_image[[si]]
      if (!length(rows)) next
      ci <- coords[rows, , drop = FALSE]
      af[rows] <- reference_area_fraction(ci, rad,
                                          min(ci[, 1]), max(ci[, 1]),
                                          min(ci[, 2]), max(ci[, 2]))
    }
  }
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
    w_glb <- exp(-d_loc / h_tech) / af[i_glb]
    ii_all[[si]] <- i_glb
    jj_all[[si]] <- j_glb
    ww_all[[si]] <- w_glb
  }
  ii <- unlist(ii_all, use.names = FALSE)
  jj <- unlist(jj_all, use.names = FALSE)
  ww <- unlist(ww_all, use.names = FALSE)
  W <- Matrix::sparseMatrix(i = ii, j = jj, x = ww, dims = c(n, n))
  W <- methods::as(W, "CsparseMatrix")
  list(W = W, image_idx = as.integer(image), n_images = nlevels(image))
}

# Detection rate per focal type, as .compute_data_informed_weights() computed it.
reference_detection_rate <- function(Y, celltype, focals) {
  ct_chr <- as.character(celltype)
  det_rate <- matrix(0, length(focals), ncol(Y), dimnames = list(focals, colnames(Y)))
  for (f in focals) {
    cells <- which(ct_chr == f)
    if (!length(cells)) next
    det_rate[f, ] <- colMeans(Y[cells, , drop = FALSE] > 0)
  }
  det_rate
}
