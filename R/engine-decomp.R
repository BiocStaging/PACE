## pace_single_frame_decomp.R — CANONICAL single-frame per-(gene,focal) decomposition.
##
## One additive frame (sums to 100%) holding cell type AND the within-cell-type blocks together,
## following Elijah's formula. Magnitudes come from OBSERVED expression (log1p CP10k); the within-part
## is split by the FIT's locked within-cell-type proportions (no re-fit lm -> respects the
## "PACE spatial from the fit" rule, [[feedback_pace_model_outputs_only_2026-05-31]]).
##
## For focal cell type c and gene g, with y_i = log1p(CP10k) observed expression:
##   denominator  = sum_{i in c} (y_i - global_mean_g)^2                      [observed]
##   Cell type    = n_c * (focal_mean_{c,g} - global_mean_g)^2               [observed lineage shift]
##   within-total = sum_{i in c} (y_i - focal_mean_{c,g})^2                  [observed]
##     Spatial    = within-total * V_state / (V_state+V_spill+V_disp)        [fit proportions]
##     Spillover  = within-total * V_spill / (...)                           [fit proportions]
##     Residual   = within-total * V_disp  / (...)                           [fit proportions]
## global_mean_g = mean over ALL cells of y; focal_mean = mean over focal cells.
##
## Cell type is the observed group-mean deviation from the tissue mean (a between-vs-global SS); the
## within blocks reuse the locked gene_focal_4block decomposition. Cell-type % is scale-dependent
## (computed on log1p CP10k; raw counts and link scale differ).

## The observed sums come from one pass over the sparse counts in C++
## (pace_single_frame_statistics_cpp): per (type, gene) the focal mean and the
## within sum of squares, and per gene the global mean, with
##   SS_within = sum over stored entries of (y - m)^2 + (n_c - stored) * m^2,
## so no dense n x G log1p matrix is built. `threads` only splits the work.
single_frame_decomp_obs <- function(Y, celltype, nCount, gene_focal_block, threads = 1L) {
  genes <- colnames(Y)
  ct <- as.character(celltype)
  TYPES <- sort(unique(ct))
  counts <- .pace_as_dgc(Y)
  frame <- pace_single_frame_statistics_cpp(counts, as.numeric(nCount), .pace_codes(ct, TYPES),
                                            length(TYPES), n_threads = .pace_thread_count(threads))
  global_mean <- frame$global_mean
  names(global_mean) <- genes
  gf <- gene_focal_block
  # Condition cohorts split the spatial component into baseline + responder
  # (gene_focal_5block has V_state_responder); no-condition cohorts have V_state.
  has_resp <- "V_state_responder" %in% names(gf)

  out <- vector("list", length(TYPES))
  k <- 0L
  for (type_index in seq_along(TYPES)) {
    fc_type <- TYPES[type_index]
    n_fc <- sum(ct == fc_type)
    if (n_fc < 5) next
    focal_mean <- frame$focal_mean[type_index, ]
    SS_within  <- frame$within_ss[type_index, ]                  # per gene, observed
    SS_lineage <- n_fc * (focal_mean - global_mean)^2            # per gene, observed
    den <- SS_lineage + SS_within
    ## fit within-proportions (locked gene_focal block) for this focal
    g4 <- gf[gf$focal == fc_type, , drop = FALSE]
    idx <- match(genes, g4$gene)
    vstate <- if (has_resp) g4$V_state_baseline[idx] else g4$V_state[idx]
    vresp  <- if (has_resp) g4$V_state_responder[idx] else 0
    vspill <- g4$V_spill[idx]
    vdisp <- g4$V_disp[idx]
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
    row$SS_lineage <- SS_lineage
    row$SS_within <- SS_within
    row$denom <- den
    row$n_focal <- n_fc
    out[[k]] <- row
  }
  do.call(rbind, out)
}

## ---------------------------------------------------------------------------
## Per-gene, per-focal variance decomposition from per-cell-type statistics.
##
## Same table as mvpql_variance_decomposition_multi() (engine-mvpql.R), without
## reading any n x G matrix. For focal type c with cells i in c, kernel columns
## N_i (K-vector), slope BLUPs u_cg with standard errors se_cg:
##
##   V_state_baseline(c, g) = var_i(N_i' u_cg) + sum_t se_tcg^2 var_i(N_it)
##                          = u_cg' S_c u_cg + diag(S_c)' se_cg^2,     S_c = Cov_{i in c}(N_i)
##   V_state_responder(c,g) = r_cg' S^R_c r_cg + diag(S^R_c)' se_rcg^2, S^R_c = Cov_{i in c}(R_i N_i[keep])
##   V_spill(c, g)          = toff_var[c, g]                                     (per-cell contamination)
##                          = beta_near,g' Cov_{i in c}(X_near,i) beta_near,g    (legacy _near fixed effects)
##   V_disp(c, g)           = log(1 + (1 + max(alpha_g, 0)) / max(mu_mean[c, g], 1e-9))   (NB1, Leckie 2020)
##
## (var(a' x) = a' Cov(x) a for a fixed vector a.) The covariances come from C++
## with stats::cov semantics; mu_mean, toff_var and the use_bleed flag come from
## .pace_fit_statistics(); count means from the sparse counts. Kept as-is from the
## original: the NB1 V_disp formula (also for NB2 fits) and its 1e-9 floor.
## ---------------------------------------------------------------------------
mvpql_variance_decomposition_stats <- function(fit, stats, df, Y, vars, X_fixed,
                                               resp_term = NULL,
                                               focal_levels = NULL,
                                               disp_model = c("nb1", "nb2"),
                                               weight_by_spec_sq = TRUE,
                                               threads = 1L) {
  disp_model <- match.arg(disp_model)
  re <- fit$re_meta
  blk_idx_ct <- which(vapply(re$blocks, `[[`, character(1), "group_col") == "celltype")
  if (length(blk_idx_ct) != 1L)
    stop("Expected exactly one celltype RE block")
  blk_ct <- re$blocks[[blk_idx_ct]]
  groups <- blk_ct$group_levels
  if (is.null(focal_levels)) focal_levels <- intersect(as.character(unique(df$celltype)), groups)
  gene_names <- colnames(fit$U)
  G <- length(gene_names)
  celltype <- as.character(df$celltype)

  ct_col <- function(t_idx, g_idx) {
    blk_ct$col_offset + (t_idx - 1L) * blk_ct$K_groups + g_idx
  }
  term2t <- stats::setNames(seq_along(blk_ct$term_levels), blk_ct$term_levels)
  has_resp <- !is.null(resp_term) && any(grepl(paste0("^", resp_term, ":"), blk_ct$term_levels))
  fix_names <- rownames(fit$B)
  spill_idx <- which(grepl("_near$|spill", fix_names, ignore.case = TRUE))
  use_bleed <- isTRUE(stats$toff_any_nonzero)

  ## per-type mean counts; a type with no cells keeps 0, as the original loop did
  group_sizes <- vapply(groups, function(type) sum(celltype == type), numeric(1))
  ct_means <- pace_group_column_means_cpp(.pace_as_dgc(Y), .pace_codes(celltype, groups),
                                          length(groups), detection = FALSE,
                                          n_threads = .pace_thread_count(threads))
  ct_means[group_sizes == 0, ] <- 0
  dimnames(ct_means) <- list(groups, gene_names)

  ## per-type covariances of the kernel columns, the responder products and the
  ## legacy spillover covariates
  kernel_cov <- .pace_group_covariances(df[, vars, drop = FALSE], celltype, groups, threads)
  responder_cov <- NULL
  if (has_resp) {
    if (!".resp_dummy" %in% colnames(df))
      stop("variance decomposition: resp_term given and interaction terms present, but df$.resp_dummy (the condition 0/1 indicator) is missing. The builder must set it.")
    responder_cov <- .pace_group_covariances(as.matrix(df[, vars, drop = FALSE]) * df[[".resp_dummy"]],
                                             celltype, groups, threads)
  }
  spill_cov <- NULL
  if (!use_bleed && length(spill_idx)) {
    spill_cov <- .pace_group_covariances(X_fixed[, spill_idx, drop = FALSE], celltype, groups, threads)
  }

  rows_naka <- vector("list", length(focal_levels))
  n_blocks <- 0L
  for (c_name in focal_levels) {
    c_idx <- match(c_name, groups)
    n_c <- sum(celltype == c_name)
    if (n_c < 5L) next
    fmean      <- ct_means[c_idx, ]
    other_max  <- apply(ct_means[-c_idx, , drop = FALSE], 2, max, na.rm = TRUE)
    spec_focal <- fmean / pmax(fmean + other_max, 1e-9)
    foratio    <- fmean / pmax(other_max, 1e-9)

    ## Baseline spatial state: quadratic form in the slope BLUPs.
    slope_rows <- vapply(vars, function(v) ct_col(term2t[[v]], c_idx), numeric(1))
    S_c <- kernel_cov[[c_name]]
    U_slopes <- fit$U[slope_rows, , drop = FALSE]                  # K x G
    SE2_slopes <- fit$se_U[slope_rows, , drop = FALSE]^2
    V_state_baseline <- colSums(U_slopes * (S_c %*% U_slopes)) +
                        colSums(SE2_slopes * diag(S_c), na.rm = TRUE)

    ## Responder spatial state: the same form over the responder interaction terms.
    V_state_responder <- rep(0, G)
    if (has_resp) {
      resp_term_names <- paste0(resp_term, ":", vars)
      matched <- vapply(resp_term_names, function(nm) {
        v <- term2t[[nm]]
        if (is.null(v)) NA_integer_ else as.integer(v)
      }, integer(1))
      keep_v <- which(!is.na(matched))
      if (length(keep_v)) {
        resp_rows <- vapply(matched[keep_v], function(ti) ct_col(ti, c_idx), numeric(1))
        S_r <- responder_cov[[c_name]][keep_v, keep_v, drop = FALSE]
        R_slopes <- fit$U[resp_rows, , drop = FALSE]
        RSE2 <- fit$se_U[resp_rows, , drop = FALSE]^2
        V_state_responder <- colSums(R_slopes * (S_r %*% R_slopes)) +
                             colSums(RSE2 * diag(S_r), na.rm = TRUE)
      }
    }

    ## Spillover: per-cell contamination offset variance, else legacy _near effects.
    V_spill <- if (use_bleed) {
      as.numeric(stats$toff_var[c_name, gene_names])
    } else if (length(spill_idx)) {
      beta_spill <- fit$B[spill_idx, , drop = FALSE]
      colSums(beta_spill * (spill_cov[[c_name]] %*% beta_spill))
    } else rep(0, G)

    mu_bar_c <- as.numeric(stats$mu_mean[c_name, gene_names])
    alpha_g <- unname(fit$alpha)
    V_disp <- if (disp_model == "nb1") {
      log(1 + (1 + pmax(alpha_g, 0)) / pmax(mu_bar_c, 1e-9))
    } else {
      log(1 + 1 / pmax(mu_bar_c, 1e-9) + pmax(alpha_g, 0))
    }

    intercept_row <- ct_col(term2t[["(Intercept)"]], c_idx)
    ct_offset_sq <- fit$U[intercept_row, ]^2 + fit$se_U[intercept_row, ]^2

    tot_5 <- ct_offset_sq + V_state_baseline + V_state_responder + V_spill + V_disp
    keep <- is.finite(tot_5) & tot_5 > 0
    spec_g <- ifelse(is.finite(spec_focal), spec_focal, 0)
    n_blocks <- n_blocks + 1L
    rows_naka[[n_blocks]] <- tibble::tibble(
      gene = gene_names, focal = c_name, n_focal = n_c,
      V_state_baseline  = unname(V_state_baseline),
      V_state_responder = unname(V_state_responder),
      V_spill = unname(V_spill), V_disp = unname(V_disp),
      celltype_offset_sq = unname(ct_offset_sq),
      focal_mean = unname(fmean), max_other_mean = unname(other_max),
      spec = unname(spec_g), focal_other_ratio = unname(foratio),
      is_contaminated = unname(!is.finite(foratio) | foratio < 1),
      `Cell type %`               = unname(ct_offset_sq      / tot_5 * 100),
      `Spatial state %`           = unname(V_state_baseline  / tot_5 * 100),
      `Responder spatial state %` = unname(V_state_responder / tot_5 * 100),
      `Spillover %`               = unname(V_spill           / tot_5 * 100),
      `Residual %`                = unname(V_disp            / tot_5 * 100)
    )[keep, ]
  }
  gene_focal_5block <- dplyr::bind_rows(rows_naka[seq_len(n_blocks)])
  .pace_decomposition_aggregates(gene_focal_5block, weight_by_spec_sq, disp_model)
}

## The aggregate tables of the decomposition (unchanged definitions from
## mvpql_variance_decomposition_multi()).
.pace_decomposition_aggregates <- function(gene_focal_5block, weight_by_spec_sq, disp_model) {
  agg_focal_5block_mean <- gene_focal_5block |>
    dplyr::group_by(.data$focal) |>
    dplyr::summarise(
      `Cell type %`               = mean(.data[["Cell type %"]],               na.rm = TRUE),
      `Spatial state %`           = mean(.data[["Spatial state %"]],           na.rm = TRUE),
      `Responder spatial state %` = mean(.data[["Responder spatial state %"]], na.rm = TRUE),
      `Spillover %`               = mean(.data[["Spillover %"]],               na.rm = TRUE),
      `Residual %`                = mean(.data[["Residual %"]],                na.rm = TRUE),
      n_genes = dplyr::n(), .groups = "drop"
    )

  agg_focal_5block_specw <- if (isTRUE(weight_by_spec_sq)) {
    gene_focal_5block |>
      dplyr::group_by(.data$focal) |>
      dplyr::summarise(
        spec_w_sum = sum(.data$spec^2, na.rm = TRUE),
        `Cell type %`               = sum(.data[["Cell type %"]]               * .data$spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        `Spatial state %`           = sum(.data[["Spatial state %"]]           * .data$spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        `Responder spatial state %` = sum(.data[["Responder spatial state %"]] * .data$spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        `Spillover %`               = sum(.data[["Spillover %"]]               * .data$spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        `Residual %`                = sum(.data[["Residual %"]]                * .data$spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        n_genes = dplyr::n(),
        n_specific_genes = sum(!.data$is_contaminated, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::select(-"spec_w_sum")
  } else NULL

  gene_focal_5block <- gene_focal_5block |>
    dplyr::mutate(Total =
      .data$V_state_baseline + .data$V_state_responder + .data$V_spill + .data$V_disp +
      .data$celltype_offset_sq)

  agg_focal_5block_pooled <- gene_focal_5block |>
    dplyr::group_by(.data$focal) |>
    dplyr::summarise(
      total_SS = sum(.data$Total, na.rm = TRUE),
      `Cell type %`               = 100 * sum(.data$celltype_offset_sq, na.rm = TRUE) / pmax(sum(.data$Total, na.rm = TRUE), 1e-12),
      `Spatial state %`           = 100 * sum(.data$V_state_baseline,   na.rm = TRUE) / pmax(sum(.data$Total, na.rm = TRUE), 1e-12),
      `Responder spatial state %` = 100 * sum(.data$V_state_responder,  na.rm = TRUE) / pmax(sum(.data$Total, na.rm = TRUE), 1e-12),
      `Spillover %`               = 100 * sum(.data$V_spill,            na.rm = TRUE) / pmax(sum(.data$Total, na.rm = TRUE), 1e-12),
      `Residual %`                = 100 * sum(.data$V_disp,             na.rm = TRUE) / pmax(sum(.data$Total, na.rm = TRUE), 1e-12),
      n_genes = dplyr::n(),
      .groups = "drop"
    )

  list(
    gene_focal_5block       = gene_focal_5block,
    agg_focal_5block_mean   = agg_focal_5block_mean,
    agg_focal_5block_specw  = agg_focal_5block_specw,
    agg_focal_5block_pooled = agg_focal_5block_pooled,
    disp_model              = disp_model
  )
}
