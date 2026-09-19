## pace_mvpql_streaming.R -- MEMORY-BOUNDED ("streaming") port of
## fit_pace_mvpql_joint_multi (the dense per-cell-HC contamination solver).
##
## GOAL: numerically identical results to fit_pace_mvpql_joint_multi for the
## canonical per-cell-HC / E^tech path, but with peak memory O(n_cells x
## chunk_size) instead of O(n_cells x n_genes).
##
## KEY IDENTITY that makes this exact:
##   The dense builder forms a dense n x G ambient matrix
##     E_tech[i, g] = (1/edge_frac[i]) *
##                    sum_{j in frNN(i, 3*h_tech), celltype(j)!=celltype(i)}
##                        exp(-d_ij / h_tech) * Y[j, g].
##   Define a SPARSE n x n weight matrix W with
##     W[i, j] = exp(-d_ij / h_tech) / edge_frac[i]
##   for cross-celltype neighbours j within 3*h_tech (0 otherwise; self
##   excluded). Then EXACTLY  E_tech[, g] = W %*% Y[, g], and for any gene
##   chunk  a_chk = W %*% Y[, chunk]  (sparse W times sparse Y-columns ->
##   dense n x |chunk|). So the full dense E_tech is NEVER materialised; we
##   stream a_chk per chunk. W is tiny because the E^tech radius is only
##   3*h_tech (a handful of neighbours per cell).
##
## SCOPE of this port (everything else identical math to the dense solver):
##   - IRLS inner solve only (USE_LAPLACE = FALSE assumed; no TMB path).
##   - additive_active (per-cell contamination) path is the supported path.
##   - No per-gene RE block (the canonical BC / Mel per-cell-HC fits have none).
##     If a per-gene block is detected we stop loudly rather than silently
##     diverge from the dense solver.
##
## All numerics are delegated to the existing helpers (build_random_design_multi,
## .solve_genes_chunk_multiblock, .alpha_nb1_mle/.alpha_nb2_mle); nothing is
## reimplemented.

fit_pace_mvpql_streaming <- function(Y, X_fixed, df, re_specs,
                                     ## the random-effect design, when the caller
                                     ## has already built it from `re_specs`
                                     re               = NULL,
                                     offset_vec       = NULL,
                                     n_iter           = 16, tol = 5e-3,
                                     ## ---- Observation family ----
                                     ## "nb"       : negative binomial with a log link, the canonical
                                     ##              count path, bit-for-bit unchanged from before.
                                     ## "gaussian" : Gaussian with an IDENTITY link, for continuous
                                     ##              intensities such as IMC protein. Working response
                                     ##              z = y and per-gene constant weight w = 1/sigma2_g.
                                     ##              Always called with contamination "none": no
                                     ##              ambient field, no rho, no anchors.
                                     family           = c("nb", "gaussian"),
                                     disp_model       = c("nb1", "nb2"),
                                     tau_shrinkage    = c("hierarchical", "shared",
                                                          "adaptive", "half_cauchy"),
                                     BPPARAM          = BiocParallel::SerialParam(),
                                     chunk_size       = 128L,
                                     alpha_max_n      = Inf,
                                     sample_weight    = NULL,
                                     n_threads        = NULL,
                                     interior_precision = 1L,
## Upper bound on the variance components. Every tau write is a pmax
                                     ## against a floor; without a ceiling a column identified only by
                                     ## the ridge climbs by mean(u^2) per iteration. 100 sits ~6-30x above
                                     ## the largest component seen on healthy cohorts (BC 12.8-17.3,
                                     ## Mel 3.2) and just below var(z) ~ 109, the variance of the working
                                     ## response the components decompose. A binding cap is reported: it
                                     ## means the term is not identified by the data.
                                     tau_max          = 100,
                                     data_informed_W  = NULL,
                                     ## ---- Streaming ambient carrier (REPLACES dense ambient_mat) ----
                                     ## ambient_W: sparse n x n dgCMatrix of E^tech weights (see header).
                                     ##   a_chk = ambient_W %*% Y[, chunk] reproduces dense E_tech[, chunk].
                                     ambient_W         = NULL,
                                     ambient_image_idx = NULL,    ## n-vector, 1-indexed image group per cell
                                     ambient_n_images  = 0L,
                                     bleed_percell     = FALSE,
                                     percell_anchor_mask = NULL,  ## n_types x G (with percell_anchor_idx) OR n x G
                                     percell_anchor_idx  = NULL,  ## length-n celltype index 1..n_types
                                     ## DEPRECATED. The fit stores the per-cell-type
                                     ## statistics every readout needs, and mu is a
                                     ## deterministic function of what it already holds
                                     ## (.pace_mu_block()), so keeping the n x G matrices
                                     ## only costs memory. TRUE still returns them, with a
                                     ## warning, for code that has not moved over yet.
                                     return_mu         = FALSE,
                                     ## ---- Guarded speed approximations (defaults ENABLE the safe mode) ----
                                     ## alpha_warmup: only re-fit the alpha dispersion MLE for the first
                                     ##   `alpha_warmup` iterations (and the last iteration); afterwards
                                     ##   freeze alpha at its warmed-up value. The alpha MLE is ~37% of
                                     ##   per-iter cost (serial) and converges quickly. Set Inf to always
                                     ##   update (exact).
                                     alpha_warmup      = 10L,
                                     ## alpha_zero_collapse: evaluate the NB1 dispersion
                                     ##   likelihood with the zero-count block collapsed into
                                     ##   one term. log f(0; mu/a, mu) = -(mu/a) log1p(a) depends
                                     ##   on the cell only through mu, so the whole zero block is
                                     ##   one term -- exact algebra, but the reassociated sum moves
                                     ##   alpha in its last digits (~1e-7 relative, measured), and
                                     ##   the fit with it. On a sparse panel ~85% of a gene's cells
                                     ##   are zeros, which is why it is worth doing. Gated on both
                                     ##   manuscript cohorts at full scale: breast cancer 1,637 and
                                     ##   melanoma 46 calls, none gained or lost, no sign flips,
                                     ##   SPP1 unchanged; fits 1.20x and 1.12x faster.
                                     alpha_zero_collapse = TRUE,
                                     ## alpha_fast_density: minimise a closed-form NB1
                                     ## objective instead of calling Rf_dnbinom_mu once
                                     ## per cell per Brent step. 5x on the dispersion
                                     ## MLE, 1.43x on a full BC fit. Agrees with the
                                     ## density it replaces to ~5e-8 on alpha; the calls
                                     ## do not move on either cohort.
                                     alpha_fast_density = TRUE,
                                     ## early_stop_tol / min_iter: break the IRLS loop once the streamed
                                     ##   MEAN rel_delta (mean over cell-genes of |Delta eta|/max(|eta|,1e-3))
                                     ##   falls below early_stop_tol, but never before min_iter iterations.
                                     ##   The mean is robust to the L-inf max, which a handful of jittering
                                     ##   cells dominate at large n (the max can oscillate forever while the
                                     ##   bulk is converged). The L-inf max is still recorded as a diagnostic.
                                     ##   Set early_stop_tol = 0 to disable.
                                     early_stop_tol    = 2e-2,
                                     min_iter          = 12L,
                                     verbose           = TRUE) {
  tau_shrinkage <- match.arg(tau_shrinkage)
  family        <- match.arg(family)
  is_gaussian   <- identical(family, "gaussian")   ## guards every Gaussian branch
  disp_model    <- match.arg(disp_model)
  if (isTRUE(return_mu))
    warning("`return_mu` is deprecated: the fit stores the per-cell-type statistics the ",
            "readouts need, and mu can be rebuilt from the fit and the counts.",
            call. = FALSE)

  ## ---- Y stays SPARSE; never densify the whole matrix. -------------------
  if (!methods::is(Y, "dgCMatrix")) {
    Y <- methods::as(methods::as(methods::as(Y, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }
  n   <- nrow(Y); g_n <- ncol(Y); p <- ncol(X_fixed)
  if (is.null(offset_vec)) offset_vec <- rep(0, n)

  ## ---- This streaming port supports IRLS only (no TMB Laplace). ----------
  INNER_SOLVE <- tolower(Sys.getenv("R_INNER_SOLVE", "irls"))
  if (!identical(INNER_SOLVE, "irls"))
    stop("fit_pace_mvpql_streaming supports only R_INNER_SOLVE=irls (got '",
         INNER_SOLVE, "'). The Laplace/polish paths are not ported.")

  if (is.null(n_threads))
    n_threads <- tryCatch(max(1L, BiocParallel::bpworkers(BPPARAM)),
                          error = function(e) 1L)
  n_threads <- max(1L, as.integer(n_threads))

  if (is.null(re)) re <- build_random_design_multi(df, re_specs)
  Z  <- re$Z; q <- ncol(Z)
  if (verbose) {
    blk_str <- paste(vapply(re$blocks, function(b)
      sprintf("%s[%dx%d]", b$group_col, b$K_terms, b$K_groups),
      character(1)), collapse = "+")
    cat(sprintf("  [mvpql.streaming] n=%d  g=%d  p_fixed=%d  q_random=%d (= %s)\n",
                n, g_n, p, q, blk_str))
  }

  ## ---- No per-gene RE block supported in the streaming port. -------------
  has_per_gene_block <- any(vapply(re$blocks,
                                   function(b) isTRUE(b$is_per_gene), logical(1)))
  if (has_per_gene_block)
    stop("fit_pace_mvpql_streaming: per-gene RE blocks are not supported by the ",
         "streaming port (canonical per-cell-HC fits have none).")

  ## ============================================================
  ## Per-cell contamination model (identifiable form):
  ## mu_ig = mu_bio_ig + rho_i * a_ig ; a_ig = E^tech ambient = (ambient_W %*% Y)[i,g].
  ## ============================================================
  percell_mode    <- isTRUE(bleed_percell)
  ## Gaussian never activates the contamination path: there is no ambient field,
  ## no rho and no anchors, because the intensities reaching it have already been
  ## corrected. add_rho stays zero throughout so the shared summary and
  ## reconstruction code, which references it, is a no-op rather than a branch.
  additive_active <- (!is_gaussian) && percell_mode &&
                     !is.null(ambient_W) && ambient_n_images > 0L
  if (!is_gaussian && !additive_active)
    stop("fit_pace_mvpql_streaming: only the per-cell contamination (additive) ",
         "path is ported. Supply bleed_percell=TRUE, ambient_W, ambient_n_images>0.")

  if (additive_active) {
    stopifnot(methods::is(ambient_W, "Matrix"),
              nrow(ambient_W) == n, ncol(ambient_W) == n,
              length(ambient_image_idx) == n)
    add_cell_image <- as.integer(ambient_image_idx)
  }
  add_rho <- numeric(n)   ## per-cell contamination loading rho_i; stays 0 on Gaussian

  if (additive_active && !is.null(percell_anchor_mask)) {
    if (!is.null(percell_anchor_idx)) {     ## MEM: per-type mask (n_types x G) + length-n index
      stopifnot(ncol(percell_anchor_mask) == g_n, length(percell_anchor_idx) == n,
                max(percell_anchor_idx) <= nrow(percell_anchor_mask))
      percell_anchor_idx <- as.integer(percell_anchor_idx)
    } else {
      stopifnot(nrow(percell_anchor_mask) == n, ncol(percell_anchor_mask) == g_n)
    }
    storage.mode(percell_anchor_mask) <- "double"
    if (verbose) cat(sprintf("    [percell_bleed] anchor mask: %s, %.0f anchor entries\n",
                             if (!is.null(percell_anchor_idx))
                               sprintf("%d types x %d genes (per-type)", nrow(percell_anchor_mask), g_n)
                             else "n x G dense",
                             sum(percell_anchor_mask)))
  }
  if (verbose)
    cat(sprintf("  [percell_bleed] activated: %d images, %d genes; mu = mu_bio + rho_i * a_ig (streaming a_ig)\n",
                ambient_n_images, g_n))

  ## ---- dispersion family ----
  disp_env <- Sys.getenv("R_DISP_MODEL", unset = "")
  disp_nb2 <- if (nzchar(disp_env)) identical(disp_env, "nb2") else identical(disp_model, "nb2")
  if (disp_nb2 && verbose) cat("  [mvpql.streaming] dispersion family = NB2 (Var=mu(1+alpha*mu))\n")

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
  ## Small, full-size state ONLY: B (p x G), U (q x G), se_B, se_U, re_var
  ## (q x G), alpha (G), tau arrays, add_rho (n). NO full n x G mu / mu_bio /
  ## mu_spill / technical_offset_mat across iterations.
  ## On the Gaussian path `alpha` carries the per-gene residual variance
  ## sigma2_g, the identity-link analogue of the NB dispersion. Seed it from the
  ## marginal per-gene variance of Y (the fit explains some of that, so it is an
  ## upper start) and let the residual-variance pass take over from iteration 1.
  if (is_gaussian) {
    alpha <- pace_marginal_variance_cpp(Y, 1e-8, n, g_n, .pace_thread_count(n_threads))
  } else {
    alpha <- rep(1, g_n)
  }
  prev_alpha <- alpha
  B    <- matrix(0, p, g_n); U <- matrix(0, q, g_n)
  re_var <- matrix(0, q, g_n)
  se_B <- matrix(NA_real_, p, g_n); se_U <- matrix(NA_real_, q, g_n)

  ## Warm-start eta from a count floor (dense uses mu <- pmax(Y, 0.5)). We need
  ## prev_eta only for the rel_delta convergence diagnostic; compute it on the
  ## fly per chunk during iter 1 (see below). Initialise prev_eta lazily.
  prev_eta_set <- FALSE
  prev_eta     <- NULL   ## becomes an n x G dense matrix ONLY if return_mu (else streamed)
  ## For the streaming rel_delta we accumulate a scalar max over chunks instead
  ## of holding the full eta matrix. We keep the previous iteration's eta as a
  ## small per-chunk-recomputable quantity: store prev_B / prev_U and recompute
  ## eta_chk = X B + Z U per chunk for both current and previous coefficients.
  prev_B <- NULL; prev_U <- NULL

  hist <- list(tau_blocks = list(), alpha = list(),
               rel_delta = numeric(),        ## L-inf max (diagnostic only)
               rel_delta_mean = numeric(),   ## mean over cell-genes (drives early stop)
               n_nan_genes = integer(),      ## genes left non-finite by the solver
               n_nonfinite = numeric(),     ## cell-genes hidden by the finite mask
               tau_max_seen = numeric())    ## largest variance component per iteration
  converged <- FALSE

  ## The EM update of the variance components is pace::tau_em_update: per
  ## random-effect column, max(mean over genes of u^2 + V, 1e-6). R only lays the
  ## columns back out as one matrix per block.
  .em_tau_blocks <- function(U_, V_) {
    tau_vec <- pace_tau_em_update_cpp(U_, V_)
    lapply(re$blocks, function(blk_i) {
      rng <- (blk_i$col_offset + 1L):(blk_i$col_offset + blk_i$n_cols)
      matrix(tau_vec[rng], blk_i$K_terms, blk_i$K_groups, byrow = TRUE,
             dimnames = list(blk_i$term_levels, blk_i$group_levels))
    })
  }

  ## ---- SPEED 1: cache a = ambient_W %*% Y ONCE (algebraically a_cache IS the
  ## full ambient E^tech matrix). The dense E_tech is NEVER materialised on disk
  ## or held twice; for BC the sparse W times a small dense Y yields a modest
  ## dense matrix; for WTA the product stays sparse. All call sites below become
  ## cheap cached column subsets instead of re-running the sparse multiply each
  ## of the (multiple) passes per iteration.
  ## The cache stays SPARSE: every pass below reads it one gene at a time in
  ## C++, so no dense n x |chunk| ambient block is ever formed.
  ## Gaussian has no ambient field. The core still receives an ambient argument
  ## at every call site, so hand it an all-zero sparse matrix of the right shape
  ## rather than branching each site: the Gaussian working response never reads
  ## it, and the zero matrix costs nothing to carry.
  a_cache <- if (is_gaussian) methods::new("dgCMatrix", Dim = c(as.integer(n), as.integer(g_n)),
                                           p = integer(g_n + 1L), i = integer(0), x = numeric(0))
             else .pace_as_dgc(ambient_W %*% Y)

  ## The fixed-effect contribution X_fixed %*% coef[, chunk], as the core builds
  ## it. When p == 1 (intercept-only under E^tech) the matrix multiply becomes a
  ## row-broadcast -- coef[1, gene] across cells, each row scaled by X_fixed[i, 1]
  ## -- which is algebraically identical for ANY p == 1, not only X_fixed == 1,
  ## and is why x1 is pulled out here. p > 1 uses the dense design instead.
  x1 <- if (p == 1L) as.numeric(X_fixed[, 1L]) else NULL
  x1_is_unit <- !is.null(x1) && all(x1 == 1)

  ## The anchor mask as the core reads it: either n_types x G with a per-cell
  ## type index, or n x G with none. An empty matrix means no masking.
  mask_matrix <- if (is.null(percell_anchor_mask)) matrix(0, 0, 0) else percell_anchor_mask
  mask_index  <- if (is.null(percell_anchor_idx)) integer(0) else as.integer(percell_anchor_idx)
  empty_matrix <- matrix(0, 0, 0)
  ## the diagnostic pass wants only the per-cell row sums, so every cell is in
  ## the same nominal group
  diagnostic_group <- integer(n)
  sample_weight_vec <- if (is.null(sample_weight)) numeric(0) else as.numeric(sample_weight)
  ## Genes inside a logical chunk are processed a few at a time, so the dense
  ## n x |sub-block| linear predictor is the only working matrix besides the
  ## chunk's z, w and ridge, which the chunk solve needs whole. The float/double
  ## gate still reads the whole logical chunk (max over its genes of colSums(w)).
  sub_genes <- max(1L, min(as.integer(chunk_size), 16L))
  ## eta for a gene chunk. The core writes the result once, into its own return
  ## buffer, in parallel over genes; the R form below it allocated three dense
  ## n x chunk matrices per call (the sparse product, the as.matrix() copy and
  ## the sum), which at 1.2M cells and a chunk of 64 was 626 MB each.
  ##
  ## The core sums the Z contributions in Z's column order and adds the fixed
  ## part afterwards, which is the association R used. That is bit-identical for
  ## p == 1; for p > 1 R sends X_fixed %*% B through BLAS and the two agree to
  ## about 2e-15, five orders inside the 1e-10 the fixtures are gated at.
  ## test-linear-predictor.R pins both against the R expression directly.
  x_fixed_dense <- if (p == 1L) matrix(0, 0, 0) else as.matrix(X_fixed)
  x1_or_empty   <- if (is.null(x1) || isTRUE(x1_is_unit)) numeric(0) else x1
  .eta_block <- function(coef_B, coef_U, genes) {
    pace_eta_block_cpp(x1_or_empty, isTRUE(x1_is_unit), x_fixed_dense, p,
                       coef_B, Z, coef_U, as.integer(genes),
                       .pace_thread_count(n_threads))
  }

  ## At iteration 1 there are no coefficients yet, so mu_bio is undefined. The
  ## dense solver seeds mu <- pmax(Y, 0.5) and mu_bio <- mu, mu_spill <- 0. We
  ## reproduce that EXACTLY in the first IRLS chunk solve by seeding chunk-local
  ## mu_bio_chk = pmax(Y_chk, 0.5), mu_chk = mu_bio_chk (add_rho is 0 at iter 1
  ## so mu_spill = 0). For iter > 1 we reconstruct mu_bio_chk / mu_chk from the
  ## current coefficients + add_rho + a_chk (identical to dense mu_bio / mu).
  ## The whole outer iteration runs in the core. There is no R loop left to
  ## fall back to: the one that used to live here was reachable only under
  ## fuse_rho, and fuse_rho existed ONLY to paper over the boundary crossings
  ## the R loop forced -- it folded the rho accumulation into the solve pass so
  ## eta was not rebuilt a third time per iteration. The core reuses eta within
  ## a chunk anyway, so the option bought nothing and is gone with its loop.
  ## git history keeps it at f4308a4 if the two ever need comparing again.
  if (!is.null(data_informed_W))
    stopifnot(identical(dim(data_informed_W), dim(tau_g_array)))


  {
    ## One call for the whole fit. The core holds B, U, re_var, alpha, tau and
    ## the contamination loading across the iterations, so the full-panel
    ## matrices each pass used to hand back -- and that R then copied into its
    ## own state -- are never allocated. Verbose lines are printed from the core
    ## as they happen; warnings are collected there and raised when it returns.
    driven <- pace_irls_driver_cpp(
      x1_or_empty, isTRUE(x1_is_unit), x_fixed_dense, X_fixed, p, Z,
      re$blocks, re$X_terms_list, re$cells_by_grp_list, re$cell_grp_list,
      Y, a_cache, offset_vec, sample_weight_vec, mask_matrix, mask_index,
      if (is.null(data_informed_W)) empty_matrix else data_informed_W,
      alpha, if (is.null(colnames(Y))) character(0) else colnames(Y),
      n, q, g_n,
      as.integer(n_iter), as.integer(min_iter), as.numeric(early_stop_tol),
      as.numeric(alpha_warmup), format(alpha_warmup),
      disp_nb2, is_gaussian, isTRUE(alpha_zero_collapse), as.numeric(alpha_max_n),
      isTRUE(alpha_fast_density), as.integer(interior_precision),
      as.integer(chunk_size), as.integer(sub_genes),
      if (is.null(tau_max)) NA_real_ else as.numeric(tau_max),
      match(tau_shrinkage, c("shared", "hierarchical", "adaptive", "half_cauchy")) - 1L,
      as.numeric(Sys.getenv("R_D0_MIN", unset = "1")),
      nzchar(Sys.getenv("R_REML_TAU", unset = "")),
      nzchar(Sys.getenv("R_RD_DIAG")), verbose, .pace_thread_count(n_threads))
    B      <- driven$B
    U      <- driven$U
    se_B   <- driven$se_B
    ## Raise the core's warnings HERE, with no C++ frame live. Raising them
    ## inside the binding meant that under options(warn = 2) the first one
    ## longjmped out of it, skipping the destructors that hand the counts, the
    ## ambient field and Z back from R's precious list -- they would stay
    ## protected for the rest of the session.
    for (message_text in driven$warnings) warning(message_text, call. = FALSE)
    se_U   <- driven$se_U
    alpha  <- driven$alpha
    add_rho     <- driven$rho
    tau_g_array <- driven$tau_g_array
    rownames(tau_g_array) <- colnames(Z)
    if (!is.null(colnames(Y))) colnames(tau_g_array) <- colnames(Y)
    it        <- as.integer(driven$n_iter)
    converged <- isTRUE(driven$converged)
    ## The core reports the variance components as the flat length-q vector the
    ## blocks are laid out in; laying them back out as one named matrix per block
    ## is the same reshape .em_tau_blocks() makes, and it is BY ROW -- the block
    ## slice runs term-major with the group varying fastest.
    tau_blocks_from_vector <- function(tau_vec) {
      lapply(re$blocks, function(blk_i) {
        rng <- (blk_i$col_offset + 1L):(blk_i$col_offset + blk_i$n_cols)
        matrix(tau_vec[rng], blk_i$K_terms, blk_i$K_groups, byrow = TRUE,
               dimnames = list(blk_i$term_levels, blk_i$group_levels))
      })
    }
    tau_blocks <- tau_blocks_from_vector(driven$tau_flat)
    run_iterations <- seq_len(it)
    hist$tau_blocks <- lapply(run_iterations,
                              function(i) tau_blocks_from_vector(driven$history_tau_flat[, i]))
    hist$alpha      <- lapply(run_iterations, function(i) driven$history_alpha[, i])
    hist$rel_delta      <- driven$history_rel_delta[run_iterations]
    hist$rel_delta_mean <- driven$history_rel_delta_mean[run_iterations]
    hist$n_nan_genes    <- as.integer(driven$history_n_nan_genes[run_iterations])
    hist$n_nonfinite    <- driven$history_n_nonfinite[run_iterations]
    hist$tau_max_seen   <- driven$history_tau_max_seen[run_iterations]
    rm(driven)
  }

  rownames(B)    <- colnames(X_fixed); rownames(U)    <- colnames(Z)
  rownames(se_B) <- colnames(X_fixed); rownames(se_U) <- colnames(Z)
  if (!is.null(colnames(Y))) {
    colnames(B) <- colnames(U) <- colnames(se_B) <- colnames(se_U) <- colnames(Y)
  }

  ## ---- mu / celltype-mean summaries (always streamed; small outputs) ----
  ## mu_celltype_means: n_celltypes x G colMeans of mu over cells of each focal
  ## celltype. mu_global_mean: length-G colMeans over all cells.
  ## When return_mu=TRUE additionally assemble the full n x G mu and
  ## technical_offset_mat (BC validation / downstream decomposition only).
  ct_block_idx <- which(vapply(re$blocks, `[[`, character(1), "group_col") == "celltype")
  ct_levels    <- if (length(ct_block_idx) == 1L) re$blocks[[ct_block_idx]]$group_levels
                  else sort(unique(as.character(df$celltype)))
  ct_chr       <- as.character(df$celltype)
  ct_code      <- .pace_codes(ct_chr, ct_levels)
  n_by_ct      <- stats::setNames(tabulate(ct_code + 1L, nbins = length(ct_levels)), ct_levels)

  mu_celltype_sum <- matrix(0, length(ct_levels), g_n,
                            dimnames = list(ct_levels, colnames(Y)))
  mu_global_sum   <- numeric(g_n)
  mu_full         <- if (return_mu) matrix(0, n, g_n, dimnames = list(NULL, colnames(Y))) else NULL
  toff_full       <- if (return_mu) matrix(0, n, g_n, dimnames = list(NULL, colnames(Y))) else NULL
  ## Per-cell contamination fraction, accumulated over gene chunks so that it is
  ## available without retaining the n x G matrices (return_mu = FALSE).
  contam_spill_sum <- numeric(n)
  contam_tot_sum   <- numeric(n)
  ## Per-(cell type, gene) variance of the technical offset over the type's cells,
  ## exactly as the decomposition computes it from the full matrix
  ## (stats::var(technical_offset_mat[cells, g], na.rm = TRUE)): every chunk holds
  ## all cells for its genes, so the per-type columns are complete here.
  toff_var <- matrix(NA_real_, length(ct_levels), g_n, dimnames = list(ct_levels, colnames(Y)))
  toff_any_nonzero <- FALSE

  ## Everything this chunk contributes is built and accumulated in C++
  ## (pace::final_pass_statistics): the two fitted-mean parts, the contamination
  ## log-offset, the per-type sums and offset variances, and the per-cell
  ## contamination sums. R holds only eta and the ambient chunk.
  chk_starts <- seq.int(1L, g_n, by = max(1L, as.integer(chunk_size)))
  for (cs in chk_starts) {
    gene_idx_chk <- cs:min(cs + chunk_size - 1L, g_n)
    eta_chk <- .eta_block(B, U, gene_idx_chk)
    pass <- pace_final_pass_statistics_cpp(eta_chk, offset_vec, a_cache, gene_idx_chk[1L],
                                           add_rho, ct_code, length(ct_levels),
                                           want_groups = TRUE, return_matrices = return_mu)
    mu_celltype_sum[, gene_idx_chk] <- pass$mu_group_sum
    toff_var[, gene_idx_chk] <- pass$toff_variance
    toff_any_nonzero <- toff_any_nonzero || pass$any_nonzero
    mu_global_sum[gene_idx_chk] <- pass$mu_column_sum
    contam_spill_sum <- contam_spill_sum + pass$spill_row_sum
    contam_tot_sum   <- contam_tot_sum   + pass$total_row_sum
    if (return_mu) {
      mu_full[, gene_idx_chk]   <- pass$mu
      toff_full[, gene_idx_chk] <- pass$toff
    }
    rm(eta_chk, pass)
  }
  mu_celltype_means <- mu_celltype_sum / pmax(n_by_ct, 1L)
  mu_global_mean    <- mu_global_sum / n
  ## Same definition as the [percell_bleed] fitting-trace diagnostic above.
  contam_frac       <- contam_spill_sum / pmax(contam_tot_sum, 1e-9)

  list(B = B, U = U, se_B = se_B, se_U = se_U,
       alpha          = alpha,
       tau_blocks     = tau_blocks,
       tau_g_array    = tau_g_array,
       tau_shrinkage  = tau_shrinkage,
       re_meta        = re,
       ## Full matrices only when return_mu (BC validation). Otherwise NULL to
       ## stay memory-bounded; mu summaries below are always available.
       mu                   = mu_full,
       technical_offset_mat = toff_full,
       bleed_offset_mat     = toff_full,
       mu_celltype_means    = mu_celltype_means,
       mu_global_mean       = mu_global_mean,
       ## Per-gene block outputs (always NULL here; no per-gene block supported).
       bleed_re_U = NULL, bleed_re_se_U = NULL,
       bleed_re_group_levels = NULL, bleed_re_cell_group = NULL,
       percell_bleed_rho = add_rho,
       contam_frac       = contam_frac,
       ## per-cell-type statistics read by paceDecompose() and paceDrivers()
       stats = .pace_statistics_record(ct_levels, n_by_ct, mu_celltype_means,
                                       toff_var, toff_any_nonzero),
       n_iter = it, converged = converged, history = hist)
}
