## pace_mvpql_joint_multi.R -- multi-block multivariate PQL fit with the
## per-cell contamination model (the identifiable bleed form).
##
## Contamination is modelled as mu_ig = mu_bio_ig + rho_i * a_ig, where a_ig is
## the RAW E^tech local cross-type ambient (ambient_mat) and rho_i is ONE
## shrunken cell-level loading shared across genes (gene-direction fixed to the
## ambient). Enabled via bleed_percell = TRUE. The legacy Φ-EM joint-bleed
## machinery (and the superseded additive-count s_g path) were removed
## 2026-06-15; both were unreachable in the canonical per-cell-HC path.

fit_pace_mvpql_joint_multi <- function(Y, X_fixed, df, re_specs,
                                       offset_vec       = NULL,
                                       pre_offset_mat   = NULL,        ## n × G PRE-COMPUTED per-cell-per-gene log-additive offset, e.g. expression-weighted technical-bleed term β_spill,g × E^tech[i,g]. Initialised into technical_offset_mat (used as a static offset; not updated unless a per-gene RE block is active).
                                       n_iter           = 16, tol = 5e-3,
                                       disp_model       = c("nb1","nb2"),  ## NB dispersion family: NB1 Var=mu(1+alpha) (default) or NB2 Var=mu(1+alpha*mu). Env R_DISP_MODEL overrides if set.
                                       tau_shrinkage    = c("hierarchical","shared","adaptive","half_cauchy"),
                                       BPPARAM          = BiocParallel::SerialParam(),
                                       chunk_size       = 128L,
                                       alpha_max_n      = Inf,
                                       sample_weight    = NULL,
                                       n_threads        = NULL,
                                       interior_precision = 1L,
                                       ## Upper bound on the variance components; a binding cap means the
                                       ## term is identified only by the ridge. See fit_pace_mvpql_streaming.
                                       tau_max          = 100,
                                       data_informed_W  = NULL,        ## OPTIONAL q × G multiplicative weight matrix for tau_g_array.  When supplied (via R_DATA_INFORMED_TAU=1 in the builder), shrinks the prior variance for (focal, neighbour, gene) random-slope columns whose data support is weak (low detection rate × low K-variance).  See .compute_data_informed_weights().
                                       ## Ambient field carriers for the per-cell contamination model (bleed_percell).
                                       ## ambient_mat (n × G) = RAW E^tech local cross-type ambient a_ig; ambient_image_idx
                                       ## / ambient_n_images give the per-image grouping. (Formerly the rank-1 joint-bleed
                                       ## args; rank-1 path removed 2026-06-15.)
                                       ambient_mat = NULL,
                                       ambient_image_idx  = NULL,   ## n-vector, 1-indexed image group per cell
                                       ambient_n_images   = 0L,
                                       bleed_percell          = FALSE,  ## PER-CELL contamination form (the identifiable one): mu_ig = mu_bio_ig + rho_i * a_ig, where a_ig = RAW E^tech (local cross-type ambient, in ambient_mat) and rho_i is ONE shrunken cell-level loading shared across all genes (fixed gene-direction = ambient). rho_i estimated by across-GENE regression per cell + empirical-Bayes shrinkage (avoids the Neyman-Scott free-per-cell trap). V_spill = Var_i(log1p(rho_i*a_ig/mu_bio)). Reuses the additive partial-offset IRLS.
                                       percell_anchor_mask    = NULL,   ## OPTIONAL anchor mask. Either n x G 0/1 (legacy), OR a memory-efficient n_types x G 0/1 mask paired with percell_anchor_idx (length-n celltype index): row tt holds the anchor genes for celltype tt, expanded by index inside the update. Identical result, ~n/n_types less memory.
                                       percell_anchor_idx     = NULL,   ## OPTIONAL length-n integer celltype index (1..n_types). When supplied, percell_anchor_mask is treated as n_types x G and expanded per-cell on the fly (avoids the dense n x G anchor matrix).
                                       verbose          = TRUE) {
  tau_shrinkage <- match.arg(tau_shrinkage)
  disp_model <- match.arg(disp_model)
  if (is.null(offset_vec)) offset_vec <- rep(0, nrow(Y))
  if (is.null(n_threads))
    n_threads <- tryCatch(max(1L, BiocParallel::bpworkers(BPPARAM)),
                           error = function(e) 1L)
  n_threads <- max(1L, as.integer(n_threads))

  ## ---- Inner solve choice: IRLS (default) vs TMB Laplace ----
  ## R_INNER_SOLVE in {"irls" (default), "laplace", "laplace_polish"}.
  ## - "irls"           PACE PQL only (legacy; biased on sparse counts)
  ## - "laplace"        TMB Laplace at every outer iteration (slow but most correct)
  ## - "laplace_polish" PACE PQL for outer iterations + ONE final TMB Laplace pass
  ##                    that corrects β, u to the marginal-likelihood maximum given
  ##                    PQL's converged τ, α.  ~5-10 min total Mel.  Recommended.
  INNER_SOLVE  <- tolower(Sys.getenv("R_INNER_SOLVE", "irls"))
  USE_LAPLACE  <- identical(INNER_SOLVE, "laplace")
  USE_POLISH   <- identical(INNER_SOLVE, "laplace_polish")
  if (USE_LAPLACE || USE_POLISH) {
    if (!exists("pace_laplace_chunk", mode = "function"))
      stop("R_INNER_SOLVE=", INNER_SOLVE, " requires pace_laplace_solver.R + the TMB DLL; ",
           "neither pace_laplace_chunk nor ensure_laplace_dll() is in the search path. ",
           "Confirm zzz.R sourced helpers/pace_laplace_solver.R.")
    get("ensure_laplace_dll", mode = "function")()
    if (verbose) cat(sprintf("  [mvpql.joint.multi] INNER SOLVE: %s\n",
                              toupper(INNER_SOLVE)))
  }
  Y <- as.matrix(Y); storage.mode(Y) <- "double"
  n <- nrow(Y); g_n <- ncol(Y); p <- ncol(X_fixed)

  re <- build_random_design_multi(df, re_specs)
  Z  <- re$Z; q <- ncol(Z)
  if (verbose) {
    blk_str <- paste(vapply(re$blocks, function(b)
      sprintf("%s[%dx%d]", b$group_col, b$K_terms, b$K_groups),
      character(1)), collapse = "+")
    cat(sprintf("  [mvpql.joint.multi] n=%d  g=%d  p_fixed=%d  q_random=%d (= %s)\n",
                n, g_n, p, q, blk_str))
  }

  ## (Rank-1 joint bleed path removed 2026-06-15 — superseded by per-cell-HC.
  ##  The ambient_mat/ambient_image_idx/ambient_n_images args now feed only the
  ##  per-cell contamination model below.)

  ## ============================================================
  ## Per-cell contamination model (the identifiable form):
  ## mu_ig = mu_bio_ig + rho_i * a_ig, where a_ig = RAW E^tech (local cross-type
  ## ambient, in ambient_mat) and rho_i is ONE shrunken cell-level loading shared
  ## across all genes (fixed gene-direction = ambient). Uses partial-offset IRLS.
  ## ============================================================
  percell_mode    <- isTRUE(bleed_percell)
  additive_active <- percell_mode &&
                     !is.null(ambient_mat) && ambient_n_images > 0L
  if (additive_active) {
    stopifnot(is.matrix(ambient_mat),
              nrow(ambient_mat) == n,
              ncol(ambient_mat) == g_n,
              length(ambient_image_idx) == n)
    pre_offset_mat   <- NULL
    add_E            <- ambient_mat     # RAW E^tech / local ambient a_ig (counts)
    storage.mode(add_E) <- "double"
    add_cell_image   <- as.integer(ambient_image_idx)
    add_rho          <- numeric(n)                 # per-cell contamination loading rho_i (percell mode)
    if (percell_mode && !is.null(percell_anchor_mask)) {
      if (!is.null(percell_anchor_idx)) {                          # MEM: per-type mask (n_types x G) + length-n celltype index
        stopifnot(ncol(percell_anchor_mask) == g_n, length(percell_anchor_idx) == n,
                  max(percell_anchor_idx) <= nrow(percell_anchor_mask))
        percell_anchor_idx <- as.integer(percell_anchor_idx)
      } else {
        stopifnot(nrow(percell_anchor_mask) == n, ncol(percell_anchor_mask) == g_n)
      }
      storage.mode(percell_anchor_mask) <- "double"
      if (verbose) cat(sprintf("    [percell_bleed] anchor mask: %s, %.0f anchor entries\n",
                               if (!is.null(percell_anchor_idx)) sprintf("%d types x %d genes (per-type)", nrow(percell_anchor_mask), g_n) else "n x G dense",
                               sum(percell_anchor_mask)))
    }
    if (verbose)
      cat(sprintf("  [percell_bleed] activated: %d images, %d genes; mu = mu_bio + rho_i * a_ig (per-cell shrunken loading, fixed ambient direction)\n",
                  ambient_n_images, g_n))
  }

  ## ---- Per-block tau matrix ----
  tau_blocks <- lapply(re$blocks, function(b) {
    matrix(1, b$K_terms, b$K_groups,
           dimnames = list(b$term_levels, b$group_levels))
  })
  build_tau_vec <- function() {
    out <- numeric(q)
    for (bi in seq_along(re$blocks)) {
      blk_i <- re$blocks[[bi]]
      slice <- as.numeric(t(tau_blocks[[bi]]))
      out[(blk_i$col_offset + 1L):(blk_i$col_offset + blk_i$n_cols)] <- slice
    }
    out
  }
  tau_g_array <- matrix(build_tau_vec(), nrow = q, ncol = g_n)
  rownames(tau_g_array) <- colnames(Z)
  if (!is.null(colnames(Y))) colnames(tau_g_array) <- colnames(Y)

  ## ---- IRLS state ----
  ## Dispersion family: NB1 (default, Var = mu(1+alpha)) or NB2 (Var = mu(1+alpha*mu)).
  ## Selected by the disp_model arg; env R_DISP_MODEL overrides if set (for builder env wiring).
  ## NB2 changes (a) the per-gene dispersion MLE, (b) the log-link IRLS working weight
  ## denominator (1+alpha) -> (1+alpha*mu), and (c) the per-cell count-scale precision
  ## weight. Default => NB1 byte-identical to the locked behaviour.
  disp_env <- Sys.getenv("R_DISP_MODEL", unset = "")
  disp_nb2 <- if (nzchar(disp_env)) identical(disp_env, "nb2") else identical(disp_model, "nb2")
  if (disp_nb2 && verbose) cat("  [mvpql.joint.multi] dispersion family = NB2 (Var=mu(1+alpha*mu))\n")
  mu    <- pmax(Y, 0.5)
  alpha <- rep(1, g_n)
  ## per-cell contamination: track mu_bio (biological mean) and mu_spill
  ## (count-scale spillover) separately; mu = mu_bio + mu_spill.
  if (additive_active) {
    mu_bio <- mu; mu_spill <- matrix(0, n, g_n)
  }
  prev_eta   <- log(mu) - offset_vec
  prev_alpha <- alpha
  hist <- list(tau_blocks = list(), alpha = list(), rel_delta = numeric(),
               n_nonfinite = integer(), tau_max_seen = numeric())
  converged <- FALSE

  B <- matrix(0, p, g_n); U <- matrix(0, q, g_n)
  re_var <- matrix(0, q, g_n)
  se_B <- matrix(NA_real_, p, g_n); se_U <- matrix(NA_real_, q, g_n)

  technical_offset_mat <- if (!is.null(pre_offset_mat)) {
    stopifnot(nrow(pre_offset_mat) == n,
              ncol(pre_offset_mat) == g_n)
    if (verbose)
      cat(sprintf("  [pre_offset] initialising technical_offset_mat from pre_offset_mat (range = [%.3f, %.3f])\n",
                   min(pre_offset_mat), max(pre_offset_mat)))
    as.matrix(pre_offset_mat)
  } else {
    matrix(0, n, g_n)
  }

  .em_tau_blocks <- function(U_, V_) {
    out <- vector("list", length(re$blocks))
    for (bi in seq_along(re$blocks)) {
      blk_i <- re$blocks[[bi]]
      m <- matrix(NA_real_, blk_i$K_terms, blk_i$K_groups,
                  dimnames = list(blk_i$term_levels, blk_i$group_levels))
      for (t in seq_len(blk_i$K_terms)) {
        for (c in seq_len(blk_i$K_groups)) {
          col <- blk_i$col_offset + (t - 1L) * blk_i$K_groups + c
          m[t, c] <- pmax(mean(U_[col, ]^2 + V_[col, ], na.rm = TRUE), 1e-6)
        }
      }
      out[[bi]] <- m
    }
    out
  }

  for (it in seq_len(n_iter)) {
    t_it <- Sys.time()
    last_iter <- (it == n_iter)
    lam_diag_mat <- 1 / tau_g_array

    chk_starts <- seq.int(1L, g_n, by = max(1L, as.integer(chunk_size)))
    for (cs in chk_starts) {
      gene_idx_chk <- cs:min(cs + chunk_size - 1L, g_n)
      iter_precision <- if (last_iter) 0L else as.integer(interior_precision)
      if (USE_LAPLACE) {
        ## TMB Laplace inner solve: takes RAW Y + offset, not (z, w).
        Y_chk      <- Y[, gene_idx_chk, drop = FALSE]
        bleed_chk  <- technical_offset_mat[, gene_idx_chk, drop = FALSE]
        offset_mat <- bleed_chk
        for (j in seq_len(ncol(bleed_chk))) offset_mat[, j] <- offset_mat[, j] + offset_vec
        tau_inv_chk <- 1 / pmax(tau_g_array[, gene_idx_chk, drop = FALSE], 1e-6)
        ## Warm start from previous iteration's BLUPs
        B_init <- B[, gene_idx_chk, drop = FALSE]
        U_init <- U[, gene_idx_chk, drop = FALSE]
        per_gene_chk <- get("pace_laplace_chunk", mode = "function")(
          Y_chunk         = Y_chk,
          X               = X_fixed,
          Z               = Z,
          offset_mat_chunk = offset_mat,
          tau_inv_chunk   = tau_inv_chk,
          alpha_chunk     = alpha[gene_idx_chk],
          beta_init_mat   = B_init,
          u_init_mat      = U_init,
          return_var      = TRUE,
          BPPARAM         = BPPARAM)
        for (jj in seq_along(gene_idx_chk)) {
          gi  <- gene_idx_chk[jj]
          res <- per_gene_chk[[jj]]
          B[, gi]      <- res$beta
          U[, gi]      <- res$u
          re_var[, gi] <- res$re_var
          if (last_iter) {
            se_B[, gi] <- sqrt(pmax(res$beta_var, 0))
            se_U[, gi] <- sqrt(pmax(res$re_var,   0))
          }
        }
        rm(per_gene_chk, Y_chk, bleed_chk, offset_mat, tau_inv_chk,
           B_init, U_init)
      } else {
        ## Legacy IRLS+ridge path (PACE PQL)
        if (additive_active) {
          ## Partial-offset IRLS: mu = mu_bio + mu_spill, dmu/deta = mu_bio.
          ## z = eta + (y - mu)/mu_bio ;  w = mu_bio^2 / (mu (1+alpha)).
          mu_bio_chk <- mu_bio[, gene_idx_chk, drop = FALSE]
          mu_chk     <- mu[,     gene_idx_chk, drop = FALSE]
          eta_chk    <- log(mu_bio_chk) - offset_vec
          z_chk      <- eta_chk + (Y[, gene_idx_chk, drop = FALSE] - mu_chk) / mu_bio_chk
          w_chk      <- if (disp_nb2) (mu_bio_chk^2) / (mu_chk * (1 + sweep(mu_chk, 2, alpha[gene_idx_chk], "*")))  # NB2: mu_bio^2/(mu(1+alpha*mu))
                        else          sweep((mu_bio_chk^2) / mu_chk, 2, (1 + alpha[gene_idx_chk]), "/")             # NB1: mu_bio^2/(mu(1+alpha))
          rm(eta_chk, mu_chk, mu_bio_chk)
        } else {
          mu_chk    <- mu[, gene_idx_chk, drop = FALSE]
          bleed_chk <- technical_offset_mat[, gene_idx_chk, drop = FALSE]
          eta_chk   <- log(mu_chk) - offset_vec - bleed_chk
          z_chk     <- eta_chk + (Y[, gene_idx_chk, drop = FALSE] - mu_chk) / mu_chk
          w_chk     <- if (disp_nb2) mu_chk / (1 + sweep(mu_chk, 2, alpha[gene_idx_chk], "*"))  # NB2: mu/(1+alpha*mu)
                       else          sweep(mu_chk, 2, (1 + alpha[gene_idx_chk]), "/")            # NB1: mu/(1+alpha)
          rm(eta_chk, mu_chk, bleed_chk)
        }
        if (!is.null(sample_weight)) w_chk <- w_chk * sample_weight
        lam_chk   <- lam_diag_mat[, gene_idx_chk, drop = FALSE]
        ## Keep the ridge out of single precision: below eps_float * sum(w) it
        ## is quantised away, at a threshold that falls as 1/n.
        if (iter_precision != 0L && .ridge_needs_double(lam_chk, w_chk))
          iter_precision <- 0L
        per_gene_chk <- .solve_genes_chunk_multiblock(
          X_fixed, re$X_terms_list, re$cell_grp_list, re$cells_by_grp_list,
          w_chk, z_chk, lam_chk, re$blocks,
          gene_idx = seq_along(gene_idx_chk),
          n_threads = n_threads,
          interior_precision = iter_precision,
          BPPARAM = BPPARAM)
        for (jj in seq_along(gene_idx_chk)) {
          gi  <- gene_idx_chk[jj]
          res <- per_gene_chk[[jj]]
          B[, gi]      <- res$beta
          U[, gi]      <- res$u
          re_var[, gi] <- pmax(res$Ainv_diag[p + seq_len(q)], 0)
          if (last_iter) {
            se_B[, gi] <- sqrt(pmax(res$Ainv_diag[seq_len(p)], 0))
            se_U[, gi] <- sqrt(pmax(res$Ainv_diag[p + seq_len(q)], 0))
          }
        }
        rm(per_gene_chk, z_chk, w_chk, lam_chk)
      }
    }

    eta_new <- as.matrix(X_fixed %*% B) + as.matrix(Z %*% U)

    ## (Rank-1 joint bleed update removed 2026-06-15 — superseded by per-cell-HC.)

    if (additive_active) {
      ## ----- Per-cell contamination update (the identifiable form) -----
      ## mu_bio = exp(offset + Xb + Zu); mu = mu_bio + rho_i * a_ig.
      mu_bio <- pmax(exp(eta_new + offset_vec), 1e-6)
      ## MEM: compute wcnt without the dense n x G a_mat.
      ## Count-scale NB precision weight wcnt = 1/Var: NB1 1/(mu(1+alpha)), NB2 1/(mu(1+alpha*mu)).
      mu_tot_pc <- pmax(mu_bio + mu_spill, 1e-8)
      wcnt <- if (disp_nb2) 1 / (mu_tot_pc * (1 + sweep(mu_tot_pc, 2, alpha, "*")))
              else          sweep(1 / mu_tot_pc, 2, (1 + alpha), "/")
      rm(mu_tot_pc)
      ## ----- Per-cell contamination loading (the identifiable form) -----
      ## rho_i estimated by across-GENE regression of residual on ambient a_ig,
      ## per cell, then empirical-Bayes shrunk (precision den_i toward prior den0)
      ## -> avoids the Neyman-Scott free-per-cell-parameter trap.
      ## SoupX-style: estimate rho from NEGATIVE-CONTROL (anchor) genes using RAW counts
      ## (their true expression ~0 -> all their counts are contamination = rho_i * a_ig).
      WA <- wcnt * add_E
      if (!is.null(percell_anchor_mask)) {
        if (!is.null(percell_anchor_idx)) {                      # MEM: expand per-type mask by celltype index (no dense n x G)
          for (tt in seq_len(nrow(percell_anchor_mask))) {
            rr <- which(percell_anchor_idx == tt)
            if (length(rr)) WA[rr, ] <- sweep(WA[rr, , drop = FALSE], 2, percell_anchor_mask[tt, ], "*")
          }
        } else WA <- WA * percell_anchor_mask
      }
      num   <- rowSums(WA * Y)                                   # per cell (RAW y on anchors)
      den   <- rowSums(WA * add_E)                               # per-cell precision
      rm(WA, wcnt)
      rho_raw <- pmax(ifelse(den > 1e-12, num / den, 0), 0)
      rho_bar <- sum(den * rho_raw) / pmax(sum(den), 1e-12)      # precision-weighted prior mean
      den_nz  <- den[is.finite(den) & den > 1e-12]
      if (length(den_nz)) {
        den0    <- stats::quantile(den_nz, 0.10, names = FALSE)  # LIGHT shrinkage (bottom 10% info)
        add_rho <- pmax((den * rho_raw + den0 * rho_bar) / (den + den0), 0)   # EB posterior mean
      } else {
        add_rho <- numeric(length(den))   # no anchor information (degenerate) -> no per-cell contamination
      }
      mu_spill <- pmax(sweep(add_E, 1, add_rho, "*"), 0)         # rho_i * a_ig (fixed gene-direction)
      mu <- pmax(mu_bio + mu_spill, 1e-6)
      technical_offset_mat <- log1p(mu_spill / mu_bio)
      if (verbose && it <= 3) {
        fr <- rowSums(mu_spill) / pmax(rowSums(mu), 1e-9)
        cat(sprintf("    [percell_bleed] it=%d  rho_i [%.4f,%.4f] med=%.4f  contam_frac med=%.3f q90=%.3f\n",
                    it, min(add_rho), max(add_rho), stats::median(add_rho),
                    stats::median(fr), stats::quantile(fr, 0.9)))
      }
      rm(num, den)
    } else {
      mu <- pmax(exp(eta_new + offset_vec + technical_offset_mat), 1e-6)
    }

    alpha_list <- BiocParallel::bplapply(seq_len(g_n),
                    function(gi) if (disp_nb2) .alpha_nb2_mle(Y[, gi], mu[, gi], max_n = alpha_max_n)
                                 else            .alpha_nb1_mle(Y[, gi], mu[, gi], max_n = alpha_max_n),
                    BPPARAM = BPPARAM)
    alpha <- unlist(alpha_list, use.names = FALSE)
    alpha[!is.finite(alpha)] <- prev_alpha[!is.finite(alpha)]
    alpha <- pmin(pmax(alpha, 1e-4), 50)
    prev_alpha <- alpha

    tau_blocks <- .em_tau_blocks(U, re_var)
    if (tau_shrinkage == "hierarchical") {
      for (bi in seq_along(tau_blocks)) {
        n_c_b <- vapply(re$cells_by_grp_list[[bi]], length, 0L)
        tau_blocks[[bi]] <- .shrink_tau_hierarchical(tau_blocks[[bi]], n_c_b)
      }
    }
    if (tau_shrinkage == "shared" || tau_shrinkage == "hierarchical") {
      tau_g_array[] <- rep(build_tau_vec(), times = g_n)
    } else {
      s2_mat <- U^2 + re_var
      rownames(s2_mat) <- rownames(tau_g_array)
      for (bi in seq_along(re$blocks)) {
        blk_i <- re$blocks[[bi]]
        rng <- (blk_i$col_offset + 1L):(blk_i$col_offset + blk_i$n_cols)
        s2_b <- s2_mat[rng, , drop = FALSE]
        d0_min_env <- as.numeric(Sys.getenv("R_D0_MIN", unset = "1"))
        tau_b_g <- if (tau_shrinkage == "adaptive") {
          .adaptive_tau_eb(s2_b, blk_i$K_terms, blk_i$K_groups,
                           gene_names = colnames(Y),
                           d0_min = d0_min_env,
                           p_fixed = p)        # propagate p_fixed for REML
        } else {
          .adaptive_tau_half_cauchy(s2_b, blk_i$K_terms, blk_i$K_groups,
                                     gene_names = colnames(Y))
        }
        tau_g_array[rng, ] <- tau_b_g
      }
    }

    ## Data-informed prior weighting: multiply tau_g_array by W (q × G) when
    ## supplied.  W reflects raw-data support for each (focal, neighbour, gene)
    ## random-slope coefficient; low support => small W => small tau => tight
    ## prior at zero => BLUP shrunk toward zero in the next per-gene WLS solve.
    ## This is PACE's structural fix for sign-flip artifacts on sparse-K /
    ## non-expressing-focal pairs.  See .compute_data_informed_weights().
    if (!is.null(data_informed_W)) {
      stopifnot(identical(dim(data_informed_W), dim(tau_g_array)))
      tau_g_array <- tau_g_array * data_informed_W
      ## Re-impose a numeric floor so 1/tau in the ridge stays finite
      tau_g_array <- pmax(tau_g_array, 1e-8)
    }
    ## ...and a ceiling, so a ridge-identified column cannot climb without end.
    tau_g_array <- .clamp_tau(tau_g_array, tau_max, it)
    hist$tau_max_seen[it] <- max(tau_g_array, na.rm = TRUE)

    rd_mat    <- abs(eta_new - prev_eta) / pmax(abs(prev_eta), 1e-3)
    ## na.rm below hides every gene whose linear predictor has gone non-finite,
    ## so count them and say so rather than reporting a clean metric.
    n_nonfinite <- sum(!is.finite(rd_mat))
    rel_delta   <- max(rd_mat, na.rm = TRUE)
    rm(rd_mat)
    if (n_nonfinite > 0L) {
      warning(sprintf(
        "iter %d: %d gene-cell linear predictors are non-finite and are excluded from the convergence metric.",
        it, n_nonfinite), call. = FALSE)
    }
    prev_eta <- eta_new
    hist$n_nonfinite[it]  <- n_nonfinite
    hist$tau_blocks[[it]] <- tau_blocks
    hist$alpha[[it]]      <- alpha
    hist$rel_delta[it]    <- rel_delta

    if (verbose) {
      tau_med <- vapply(tau_blocks, function(m) stats::median(m), numeric(1))
      cat(sprintf("  [mvpql.joint.multi] iter %d  rel_delta=%.3g  alpha[med]=%.2f  tau=[%s]  (%.1fs)\n",
                  it, rel_delta, median(alpha),
                  paste(sprintf("%.3f", tau_med), collapse=","),
                  as.numeric(difftime(Sys.time(), t_it, units = "secs"))))
    }
    ## Early-exit on convergence is intentionally DISABLED here: all canonical
    ## callers ran the full n_iter (the legacy guard was gated on bleed_start_iter,
    ## which was always set > n_iter, so it never fired).  Preserving that means
    ## `converged` stays FALSE and the loop always runs n_iter iterations.
    invisible(gc(verbose = FALSE))
  }

  ## ============================================================
  ## LAPLACE POLISH (R_INNER_SOLVE=laplace_polish)
  ## ============================================================
  ## After PACE PQL converges, do ONE pass over all genes with a proper TMB
  ## Laplace approximation given the converged τ, α, technical_offset_mat.  This
  ## corrects PACE's joint-mode BLUPs to the marginal-likelihood maximum that
  ## glmmTMB-style NB1 ML would produce, removing the 10-15x inflation seen
  ## on sparse multi-collinear K_<neighbour> random slopes.  τ/α/bleed_offset
  ## are NOT updated (no feedback loop) — this is a single corrective pass.
  if (USE_POLISH) {
    t_pol <- Sys.time()
    if (verbose) cat(sprintf("  [mvpql.joint.multi] Laplace polish pass over %d genes ...\n", g_n))
    ## Default to SERIAL polish: TMB autodiff tape ~500 MB per gene per worker.
    ## Each parallel worker holds its own tape, so 2 workers ≈ 1 GB extra per
    ## gene that R's lazy GC doesn't reclaim until end-of-pass.  Serial + per-
    ## gene gc() bounds peak memory.  Override with R_LAPLACE_WORKERS only if
    ## memory headroom is verified.
    wkrs <- as.integer(Sys.getenv("R_LAPLACE_WORKERS", unset = "1"))
    bp_pol <- if (wkrs > 1L)
                BiocParallel::MulticoreParam(workers = wkrs, RNGseed = 1L)
              else BiocParallel::SerialParam()
    pol_chunk_size <- max(1L, as.integer(chunk_size))
    chk_starts <- seq.int(1L, g_n, by = pol_chunk_size)
    n_done <- 0L
    for (cs in chk_starts) {
      gene_idx_chk <- cs:min(cs + pol_chunk_size - 1L, g_n)
      Y_chk        <- Y[, gene_idx_chk, drop = FALSE]
      bleed_chk    <- technical_offset_mat[, gene_idx_chk, drop = FALSE]
      off_mat_chk  <- bleed_chk
      for (j in seq_len(ncol(off_mat_chk)))
        off_mat_chk[, j] <- off_mat_chk[, j] + offset_vec
      tau_inv_chk  <- 1 / pmax(tau_g_array[, gene_idx_chk, drop = FALSE], 1e-6)
      B_init       <- B[, gene_idx_chk, drop = FALSE]
      U_init       <- U[, gene_idx_chk, drop = FALSE]
      pol_chk <- get("pace_laplace_chunk", mode = "function")(
        Y_chunk         = Y_chk,
        X               = X_fixed,
        Z               = Z,
        offset_mat_chunk = off_mat_chk,
        tau_inv_chunk   = tau_inv_chk,
        alpha_chunk     = alpha[gene_idx_chk],
        beta_init_mat   = B_init,
        u_init_mat      = U_init,
        return_var      = TRUE,
        BPPARAM         = bp_pol)
      for (jj in seq_along(gene_idx_chk)) {
        gi  <- gene_idx_chk[jj]
        res <- pol_chk[[jj]]
        if (is.null(res)) next
        B[, gi]      <- res$beta
        U[, gi]      <- res$u
        re_var[, gi] <- res$re_var
        se_B[, gi]   <- sqrt(pmax(res$beta_var, 0))
        se_U[, gi]   <- sqrt(pmax(res$re_var,   0))
      }
      rm(pol_chk, Y_chk, bleed_chk, off_mat_chk, tau_inv_chk, B_init, U_init)
      n_done <- n_done + length(gene_idx_chk)
      if (verbose) cat(sprintf("    [polish] %d / %d genes done (%.1fs elapsed)\n",
                                n_done, g_n,
                                as.numeric(difftime(Sys.time(), t_pol, units = "secs"))))
      invisible(gc(verbose = FALSE))
    }
    ## Recompute mu from polished (B, U)
    eta_pol <- as.matrix(X_fixed %*% B) + as.matrix(Z %*% U)
    mu <- pmax(exp(eta_pol + offset_vec + technical_offset_mat), 1e-6)
    rm(eta_pol)
    if (verbose) cat(sprintf("  [mvpql.joint.multi] Laplace polish done in %.1f min\n",
                              as.numeric(difftime(Sys.time(), t_pol, units = "mins"))))
  }

  rownames(B)    <- colnames(X_fixed); rownames(U)    <- colnames(Z)
  rownames(se_B) <- colnames(X_fixed); rownames(se_U) <- colnames(Z)
  if (!is.null(colnames(Y))) {
    colnames(B) <- colnames(U) <- colnames(se_B) <- colnames(se_U) <- colnames(Y)
  }

  ## Per-gene bleed RE block removed 2026-06-25; these schema fields stay NULL
  ## (kept for output compatibility with the streaming solver).
  bleed_re_U            <- NULL
  bleed_re_se_U         <- NULL
  bleed_re_group_levels <- NULL
  bleed_re_cell_group   <- NULL

  ## Per-cell contamination outputs.
  percell_bleed_rho     <- if (additive_active && percell_mode) add_rho else NULL
  mu_spill_out          <- if (additive_active) mu_spill   else NULL
  mu_bio_out            <- if (additive_active) mu_bio     else NULL
  ## Same definition as the [percell_bleed] fitting-trace diagnostic, and as the
  ## streaming solver's `contam_frac`.
  contam_frac_out       <- if (additive_active)
                             rowSums(mu_spill) / pmax(rowSums(mu), 1e-9) else NULL

  ## Per-cell-type statistics read by paceDecompose() and paceDrivers(), from the
  ## dense matrices this solver already holds.
  ct_block <- which(vapply(re$blocks, `[[`, character(1), "group_col") == "celltype")
  ct_levels <- if (length(ct_block) == 1L) re$blocks[[ct_block]]$group_levels
               else sort(unique(as.character(df$celltype)))
  stats <- .pace_statistics_from_matrices(mu, technical_offset_mat,
                                          as.character(df$celltype), ct_levels,
                                          gene_names = colnames(Y))

  list(B = B, U = U, se_B = se_B, se_U = se_U,
       alpha          = alpha,
       tau_blocks     = tau_blocks,
       tau_g_array    = tau_g_array,
       tau_shrinkage  = tau_shrinkage,
       mu             = mu,
       re_meta        = re,
       technical_offset_mat = technical_offset_mat,
       ## Backward-compat alias (kept so existing figures using `bleed_offset_mat`
       ## still work).  Both point to the same matrix.
       bleed_offset_mat     = technical_offset_mat,
       bleed_re_U            = bleed_re_U,
       bleed_re_se_U         = bleed_re_se_U,
       bleed_re_group_levels = bleed_re_group_levels,
       bleed_re_cell_group   = bleed_re_cell_group,
       percell_bleed_rho      = percell_bleed_rho,
       mu_spill               = mu_spill_out,
       mu_bio                 = mu_bio_out,
       contam_frac            = contam_frac_out,
       mu_celltype_means      = stats$mu_mean,
       stats                  = stats,
       n_iter         = it, converged = converged, history = hist)
}
