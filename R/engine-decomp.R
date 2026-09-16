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

## The observed sums come from one pass over the sparse counts
## (pace::single_frame_statistics) and the shares from pace::single_frame_shares;
## no dense n x G log1p matrix is built. R only aligns each focal's fit rows to
## the gene order and assembles the returned frame.
single_frame_decomp_obs <- function(Y, celltype, nCount, gene_focal_block, threads = 1L) {
  genes <- colnames(Y)
  n_genes <- length(genes)
  ct <- as.character(celltype)
  TYPES <- sort(unique(ct))
  threads <- .pace_thread_count(threads)
  type_code <- .pace_codes(ct, TYPES)
  group_size <- tabulate(type_code + 1L, nbins = length(TYPES))
  frame <- pace_single_frame_statistics_cpp(.pace_as_dgc(Y), as.numeric(nCount), type_code,
                                            length(TYPES), n_threads = threads)
  gf <- gene_focal_block
  # Condition cohorts split the spatial component into baseline + responder
  # (gene_focal_5block has V_state_responder); no-condition cohorts have V_state.
  has_resp <- "V_state_responder" %in% names(gf)

  ## the fit's within-cell-type components, one row per type, aligned to `genes`
  empty <- matrix(NA_real_, length(TYPES), n_genes)
  v_state <- empty
  v_responder <- if (has_resp) empty else matrix(0, 0, 0)
  v_spill <- empty
  v_disp <- empty
  for (type_index in seq_along(TYPES)) {
    g4 <- gf[gf$focal == TYPES[type_index], , drop = FALSE]
    idx <- match(genes, g4$gene)
    v_state[type_index, ] <- if (has_resp) g4$V_state_baseline[idx] else g4$V_state[idx]
    if (has_resp) v_responder[type_index, ] <- g4$V_state_responder[idx]
    v_spill[type_index, ] <- g4$V_spill[idx]
    v_disp[type_index, ] <- g4$V_disp[idx]
  }
  shares <- pace_single_frame_shares_cpp(frame$focal_mean, frame$within_ss, frame$global_mean,
                                         group_size, v_state, v_responder, v_spill, v_disp,
                                         n_threads = threads)

  ## One block per kept focal, row-named by gene, as the per-cell version was:
  ## its first numeric column carried the gene names of colMeans(), and rbind()
  ## made them unique across focals.
  kept <- which(shares$keep == 1L)
  if (!length(kept)) return(NULL)
  out <- vector("list", length(kept))
  for (k in seq_along(kept)) {
    rows <- seq_len(n_genes) + (kept[k] - 1L) * n_genes
    block <- data.frame(
      focal = TYPES[kept[k]], gene = genes,
      `Cell type %` = shares$pct_celltype[rows],
      `Spatial %`   = shares$pct_spatial[rows],
      check.names = FALSE, stringsAsFactors = FALSE, row.names = genes)
    if (has_resp) block[["Responder spatial %"]] <- shares$pct_responder[rows]
    block[["Spillover %"]] <- shares$pct_spill[rows]
    block[["Residual %"]]  <- shares$pct_residual[rows]
    block$SS_lineage <- shares$SS_lineage[rows]
    block$SS_within <- shares$SS_within[rows]
    block$denom <- shares$denom[rows]
    block$n_focal <- group_size[kept[k]]
    out[[k]] <- block
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
## (var(a' x) = a' Cov(x) a for a fixed vector a.) The arithmetic is
## pace::variance_decomposition; R locates the BLUP rows by name, hands over the
## covariances and the statistics, and assembles the table. Kept as-is from the
## original: the NB1 V_disp formula (also for NB2 fits) and its 1e-9 floor.
## ---------------------------------------------------------------------------
mvpql_variance_decomposition_stats <- function(fit, stats, df, Y, vars, X_fixed,
                                               resp_term = NULL,
                                               focal_levels = NULL,
                                               disp_model = c("nb1", "nb2"),
                                               weight_by_spec_sq = TRUE,
                                               threads = 1L) {
  disp_model <- match.arg(disp_model)
  threads <- .pace_thread_count(threads)
  re <- fit$re_meta
  blk_idx_ct <- which(vapply(re$blocks, `[[`, character(1), "group_col") == "celltype")
  if (length(blk_idx_ct) != 1L)
    stop("Expected exactly one celltype RE block")
  blk_ct <- re$blocks[[blk_idx_ct]]
  groups <- blk_ct$group_levels
  if (is.null(focal_levels)) focal_levels <- intersect(as.character(unique(df$celltype)), groups)
  gene_names <- colnames(fit$U)
  celltype <- as.character(df$celltype)

  ## 1-based row of fit$U holding term `t_idx` of group `g_idx`
  ct_col <- function(t_idx, g_idx) {
    blk_ct$col_offset + (t_idx - 1L) * blk_ct$K_groups + g_idx
  }
  term2t <- stats::setNames(seq_along(blk_ct$term_levels), blk_ct$term_levels)
  has_resp <- !is.null(resp_term) && any(grepl(paste0("^", resp_term, ":"), blk_ct$term_levels))
  fix_names <- rownames(fit$B)
  spill_idx <- which(grepl("_near$|spill", fix_names, ignore.case = TRUE))
  use_bleed <- isTRUE(stats$toff_any_nonzero)
  focal_group <- match(focal_levels, groups)
  type_code <- .pace_codes(celltype, groups)
  group_size <- tabulate(type_code + 1L, nbins = length(groups))

  ## per-type mean counts, and per-type covariances of the kernel columns, the
  ## responder products and the legacy spillover covariates
  ct_means <- pace_group_column_means_cpp(.pace_as_dgc(Y), type_code, length(groups),
                                          detection = FALSE, n_threads = threads)
  kernel_cov <- .pace_group_covariances(df[, vars, drop = FALSE], celltype, groups, threads)
  empty_matrix <- matrix(0, 0, 0)
  responder_cov <- numeric(0)
  responder_rows <- matrix(0L, 0, 0)
  responder_keep <- integer(0)
  if (has_resp) {
    if (!".resp_dummy" %in% colnames(df))
      stop("variance decomposition: resp_term given and interaction terms present, but df$.resp_dummy (the condition 0/1 indicator) is missing. The builder must set it.")
    responder_cov <- .pace_group_covariances(as.matrix(df[, vars, drop = FALSE]) * df[[".resp_dummy"]],
                                             celltype, groups, threads)
    matched <- match(paste0(resp_term, ":", vars), names(term2t))
    responder_keep <- which(!is.na(matched))
    responder_rows <- vapply(focal_group, function(g_idx)
      as.integer(ct_col(term2t[matched[responder_keep]], g_idx)), integer(length(responder_keep)))
    dim(responder_rows) <- c(length(responder_keep), length(focal_group))
  }
  spill_cov <- numeric(0)
  beta_spill <- empty_matrix
  if (!use_bleed && length(spill_idx)) {
    spill_cov <- .pace_group_covariances(X_fixed[, spill_idx, drop = FALSE], celltype, groups, threads)
    beta_spill <- fit$B[spill_idx, , drop = FALSE]
  }

  ## The statistics are indexed by position, from names resolved once, so a
  ## mismatch is an error rather than a silently repeated row.
  stat_rows <- match(focal_levels, rownames(stats$mu_mean))
  stat_cols <- match(gene_names, colnames(stats$mu_mean))
  if (anyNA(stat_rows) || anyNA(stat_cols))
    stop("the fit's per-cell-type statistics do not cover every focal cell type and gene; ",
         "re-run paceDecompose() on the object the model was fitted to.", call. = FALSE)
  statistic_block <- function(values) values[stat_rows, stat_cols, drop = FALSE]

  if (anyNA(term2t[vars]))
    stop("variance decomposition: the celltype random-effect block has no term for ",
         paste(vars[is.na(term2t[vars])], collapse = ", "), call. = FALSE)
  slope_rows <- vapply(focal_group, function(g_idx)
    as.integer(ct_col(term2t[vars], g_idx)), integer(length(vars)))
  dim(slope_rows) <- c(length(vars), length(focal_group))
  intercept_rows <- as.integer(ct_col(term2t[["(Intercept)"]], focal_group))

  blocks <- pace_variance_decomposition_cpp(
    ct_means = ct_means, group_size = group_size,
    focal_group = as.integer(focal_group), n_focal = as.integer(group_size[focal_group]),
    u = fit$U, se_u = fit$se_U, slope_rows = slope_rows, kernel_cov = kernel_cov,
    responder_rows = responder_rows, responder_keep = as.integer(responder_keep),
    responder_cov = responder_cov, intercept_rows = intercept_rows,
    toff_var = if (use_bleed) statistic_block(stats$toff_var) else empty_matrix,
    spill_cov = spill_cov, beta_spill = beta_spill,
    mu_mean = statistic_block(stats$mu_mean),
    alpha = unname(fit$alpha), nb1 = disp_model == "nb1", n_threads = threads)

  keep <- blocks$keep == 1L
  ## The per-cell implementation carried the gene names of the count means into
  ## these four columns, and older stored tables have them, so they are restored
  ## here rather than left to the vectors the core returns.
  row_genes <- rep(gene_names, times = length(focal_levels))
  named_by_gene <- function(values) stats::setNames(values[keep], row_genes[keep])
  gene_focal_5block <- tibble::tibble(
    gene = rep(gene_names, times = length(focal_levels)),
    focal = rep(focal_levels, each = length(gene_names)),
    n_focal = rep(as.integer(group_size[focal_group]), each = length(gene_names)),
    V_state_baseline = blocks$V_state_baseline, V_state_responder = blocks$V_state_responder,
    V_spill = blocks$V_spill, V_disp = blocks$V_disp,
    celltype_offset_sq = blocks$celltype_offset_sq,
    focal_mean = blocks$focal_mean, max_other_mean = blocks$max_other_mean,
    spec = blocks$spec, focal_other_ratio = blocks$focal_other_ratio,
    is_contaminated = blocks$is_contaminated == 1L,
    `Cell type %`               = blocks$pct_celltype,
    `Spatial state %`           = blocks$pct_state,
    `Responder spatial state %` = blocks$pct_responder,
    `Spillover %`               = blocks$pct_spill,
    `Residual %`                = blocks$pct_residual)[keep, ]
  gene_focal_5block$focal_mean <- named_by_gene(blocks$focal_mean)
  gene_focal_5block$max_other_mean <- named_by_gene(blocks$max_other_mean)
  gene_focal_5block$spec <- named_by_gene(blocks$spec)
  gene_focal_5block$focal_other_ratio <- named_by_gene(blocks$focal_other_ratio)
  .pace_decomposition_aggregates(gene_focal_5block, weight_by_spec_sq, disp_model)
}

## The aggregate tables of the decomposition (unchanged definitions from
## mvpql_variance_decomposition_multi()): the arithmetic is
## pace::decomposition_aggregates, and R only fixes the row order (the order
## dplyr::group_by() would have produced) and names the columns.
.pace_decomposition_aggregates <- function(gene_focal_5block, weight_by_spec_sq, disp_model) {
  block_names <- c("Cell type %", "Spatial state %", "Responder spatial state %",
                   "Spillover %", "Residual %")
  component_names <- c("celltype_offset_sq", "V_state_baseline", "V_state_responder",
                       "V_spill", "V_disp")
  focal_levels <- dplyr::group_keys(dplyr::group_by(gene_focal_5block, .data$focal))$focal
  focal_code <- match(gene_focal_5block$focal, focal_levels) - 1L
  aggregates <- pace_decomposition_aggregates_cpp(
    focal_code = as.integer(focal_code), n_focals = length(focal_levels),
    pct = as.matrix(gene_focal_5block[, block_names]),
    components = as.matrix(gene_focal_5block[, component_names]),
    spec = gene_focal_5block$spec,
    is_contaminated = as.integer(gene_focal_5block$is_contaminated),
    weight_by_spec_sq = isTRUE(weight_by_spec_sq))

  as_table <- function(values) {
    out <- tibble::tibble(focal = focal_levels)
    for (b in seq_along(block_names)) out[[block_names[b]]] <- values[, b]
    out
  }
  agg_focal_5block_mean <- as_table(aggregates$mean)
  agg_focal_5block_mean$n_genes <- aggregates$n_genes

  agg_focal_5block_specw <- NULL
  if (isTRUE(weight_by_spec_sq)) {
    agg_focal_5block_specw <- as_table(aggregates$specw)
    agg_focal_5block_specw$n_genes <- aggregates$n_genes
    agg_focal_5block_specw$n_specific_genes <- aggregates$n_specific_genes
  }

  gene_focal_5block$Total <- aggregates$Total
  agg_focal_5block_pooled <- tibble::tibble(focal = focal_levels, total_SS = aggregates$total_SS)
  for (b in seq_along(block_names))
    agg_focal_5block_pooled[[block_names[b]]] <- aggregates$pooled[, b]
  agg_focal_5block_pooled$n_genes <- aggregates$n_genes

  list(
    gene_focal_5block       = gene_focal_5block,
    agg_focal_5block_mean   = agg_focal_5block_mean,
    agg_focal_5block_specw  = agg_focal_5block_specw,
    agg_focal_5block_pooled = agg_focal_5block_pooled,
    disp_model              = disp_model
  )
}
