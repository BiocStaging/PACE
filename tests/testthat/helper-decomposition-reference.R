# The dense single-frame decomposition as implemented at f55d976 (R/engine-decomp.R),
# kept verbatim as the reference for the sparse-statistics rewrite.

reference_single_frame_decomp_obs <- function(Y, celltype, nCount, gene_focal_block) {
  genes <- colnames(Y); ct <- as.character(celltype); TYPES <- sort(unique(ct))
  Ylog <- log1p(Y * (1e4 / nCount))                      # n x G observed log1p CP10k
  global_mean <- as.numeric(Matrix::colMeans(Ylog)); names(global_mean) <- genes
  gf <- gene_focal_block
  # Condition cohorts split the spatial component into baseline + responder
  # (gene_focal_5block has V_state_responder); no-condition cohorts have V_state.
  has_resp <- "V_state_responder" %in% names(gf)

  out <- vector("list", length(TYPES)); k <- 0L
  for (fc_type in TYPES) {
    fc <- which(ct == fc_type); n_fc <- length(fc); if (n_fc < 5) next
    Yf <- as.matrix(Ylog[fc, , drop = FALSE])
    focal_mean <- colMeans(Yf)
    SS_within  <- colSums(sweep(Yf, 2, focal_mean, "-")^2)     # per gene, observed
    SS_lineage <- n_fc * (focal_mean - global_mean)^2          # per gene, observed
    den <- SS_lineage + SS_within
    ## fit within-proportions (locked gene_focal block) for this focal
    g4 <- gf[gf$focal == fc_type, , drop = FALSE]
    idx <- match(genes, g4$gene)
    vstate <- if (has_resp) g4$V_state_baseline[idx] else g4$V_state[idx]
    vresp  <- if (has_resp) g4$V_state_responder[idx] else 0
    vspill <- g4$V_spill[idx]; vdisp <- g4$V_disp[idx]
    vt <- vstate + vresp + vspill + vdisp
    p_sp   <- ifelse(vt > 0, vstate / vt, 0)
    p_resp <- ifelse(vt > 0, vresp  / vt, 0)
    p_bl   <- ifelse(vt > 0, vspill / vt, 0)
    p_rs   <- ifelse(vt > 0, vdisp  / vt, 0)
    k <- k + 1L
    row <- data.frame(
      focal = fc_type, gene = genes,
      `Cell type %` = 100 * SS_lineage       / den,
      `Spatial %`   = 100 * SS_within * p_sp / den,
      check.names = FALSE, stringsAsFactors = FALSE)
    if (has_resp) row[["Responder spatial %"]] <- 100 * SS_within * p_resp / den
    row[["Spillover %"]] <- 100 * SS_within * p_bl / den
    row[["Residual %"]]  <- 100 * SS_within * p_rs / den
    row$SS_lineage <- SS_lineage; row$SS_within <- SS_within
    row$denom <- den; row$n_focal <- n_fc
    out[[k]] <- row
  }
  do.call(rbind, out)
}

# The f55d976 decomposition wrapper: the per-cell loop of
# mvpql_variance_decomposition_multi() on a fit that stores mu and the technical
# offset, plus the 4-block view.
reference_pace_decompose <- function(fit, df, Y, types, X_fixed, resp_term = NULL) {
  dec <- PACE:::mvpql_variance_decomposition_multi(
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

# Two objects are equal when their non-numeric parts are identical and every
# numeric entry satisfies |a - b| <= tol_abs + tol_rel * |b|, with the absolute
# floor scaled to the column (exact zeros on one side, 1e-17 residue on the other).
# tol_rel = 1e-10 leaves more than 200x headroom over the measured identity errors
# (<= 4.6e-13); tol_abs = 1e-14 x max(1, max |b|) covers quadratic-form residue
# where the per-cell variance is exactly zero.
expect_tables_close <- function(new, reference, tol_rel = 1e-10, tol_abs_scale = 1e-14, label = "") {
  if (is.data.frame(reference)) {
    expect_identical(names(new), names(reference), label = label)
    expect_identical(nrow(new), nrow(reference), label = label)
    for (column in names(reference)) {
      a <- new[[column]]
      b <- reference[[column]]
      if (is.numeric(b)) {
        a <- as.numeric(a)
        b <- as.numeric(b)
        expect_identical(is.finite(a), is.finite(b), label = paste(label, column))
        finite <- is.finite(b)
        tol_abs <- tol_abs_scale * max(1, abs(b[finite]))
        within <- abs(a[finite] - b[finite]) <= tol_abs + tol_rel * abs(b[finite])
        expect_true(all(within), label = paste(label, column))
      } else {
        expect_identical(a, b, label = paste(label, column))
      }
    }
    return(invisible(TRUE))
  }
  if (is.list(reference)) {
    expect_identical(names(new), names(reference), label = label)
    for (key in names(reference)) {
      expect_tables_close(new[[key]], reference[[key]], tol_rel, tol_abs_scale, paste(label, key))
    }
    return(invisible(TRUE))
  }
  if (is.numeric(reference)) {
    expect_tables_close(data.frame(value = as.numeric(new)), data.frame(value = as.numeric(reference)),
                        tol_rel, tol_abs_scale, label)
    return(invisible(TRUE))
  }
  expect_identical(new, reference, label = label)
}
