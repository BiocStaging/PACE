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
                                     ## fuse_rho: fold the per-cell rho accumulation INTO the solve pass
                                     ##   (Pass 1), eliminating the separate second chunk pass. The rho
                                     ##   update then uses the PRE-solve (iteration t-1) B,U -- a one-step
                                     ##   lag; same fixed point (B,U are stable at convergence), different
                                     ##   path. Because the pre-solve eta equals the previous iteration's
                                     ##   eta, the eta-based rel_delta is unusable, so the fused path uses a
                                     ##   COEFFICIENT-based convergence metric (mean |dU|/max(|U|,1e-3)).
                                     ##   ~1 of 3 chunk matmuls + 1 of 2 Y densifies saved per iter.
                                     fuse_rho          = FALSE,
                                     verbose           = TRUE) {
  tau_shrinkage <- match.arg(tau_shrinkage)
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
  additive_active <- percell_mode &&
                     !is.null(ambient_W) && ambient_n_images > 0L
  if (!additive_active)
    stop("fit_pace_mvpql_streaming: only the per-cell contamination (additive) ",
         "path is ported. Supply bleed_percell=TRUE, ambient_W, ambient_n_images>0.")

  stopifnot(methods::is(ambient_W, "Matrix"),
            nrow(ambient_W) == n, ncol(ambient_W) == n,
            length(ambient_image_idx) == n)
  add_cell_image <- as.integer(ambient_image_idx)
  add_rho        <- numeric(n)   ## per-cell contamination loading rho_i

  if (!is.null(percell_anchor_mask)) {
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
  alpha      <- rep(1, g_n)
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
  a_cache <- .pace_as_dgc(ambient_W %*% Y)

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
  for (it in seq_len(n_iter)) {
    t_it      <- Sys.time()
    ## SPEED 4 (guarded): EARLY STOP. Decide at the TOP of the iteration using
    ## the PREVIOUS iteration's MEAN rel_delta so that the stopping iteration is
    ## treated as last_iter (runs SEs at precision 0 + final alpha update),
    ## keeping the output schema identical to a full run. The MEAN (not the L-inf
    ## max) is the stopping metric: at large n the max is dominated by a few
    ## jittering cells and never settles, so the max-based stop would never fire
    ## even when the bulk has converged. We never stop before min_iter.
    ## early_stop_tol = 0 disables (rel_delta_mean is always >= 0).
    ## A fit that has lost genes must not stop early and report success, so the
    ## previous iteration's non-finite gene count gates the early stop too.
    stop_now <- (it > min_iter && early_stop_tol > 0 &&
                 length(hist$rel_delta_mean) >= 1L &&
                 is.finite(hist$rel_delta_mean[it - 1L]) &&
                 hist$rel_delta_mean[it - 1L] < early_stop_tol &&
                 isTRUE(hist$n_nan_genes[it - 1L] == 0L))
    last_iter <- (it == n_iter) || stop_now
    lam_diag_mat <- 1 / tau_g_array

    ## rho + convergence accumulators (populated by Pass 1 if fuse_rho, else Pass 2)
    num <- numeric(n); den <- numeric(n)
    rel_delta <- 0; rd_n <- 0; rd_sum <- 0; rd_nonfinite <- 0
    nan_genes_iter <- integer(0)
    rd_g01 <- 0; rd_g05 <- 0; rd_g1 <- 0; rd_g10 <- 0
    RD_DIAG <- nzchar(Sys.getenv("R_RD_DIAG"))

    chk_starts <- seq.int(1L, g_n, by = max(1L, as.integer(chunk_size)))
    if (!fuse_rho) {
      ## The whole chunk loop runs in the core. The design is converted once for
      ## the pass rather than once per chunk, z and w are buffers that never
      ## become R matrices, and only the finished coefficients come back. The
      ## calls it makes, and their order, are the ones the R loop below makes.
      pass1 <- pace_fit_pass1_cpp(
        x1_or_empty, isTRUE(x1_is_unit), x_fixed_dense, X_fixed, p, Z,
        re$blocks, re$X_terms_list, re$cells_by_grp_list, re$cell_grp_list,
        B, U, lam_diag_mat, Y, a_cache, offset_vec, add_rho, alpha, sample_weight_vec,
        disp_nb2, it == 1L, n, as.integer(chunk_size), as.integer(sub_genes),
        as.integer(interior_precision), last_iter, .pace_thread_count(n_threads))
      B      <- pass1$B
      U      <- pass1$U
      re_var <- pass1$re_var
      if (last_iter) {
        se_B <- pass1$se_B
        se_U <- pass1$se_U
      }
      nan_genes_iter <- pass1$nan_genes
      if (length(nan_genes_iter) && verbose)
        cat(sprintf("    [nan-guard] it=%d: %d gene(s) STILL non-finite after double solve\n",
                    it, length(nan_genes_iter)))
      rm(pass1)
    } else {
      for (cs in chk_starts) {
        gene_idx_chk   <- cs:min(cs + chunk_size - 1L, g_n)
        iter_precision <- if (last_iter) 0L else as.integer(interior_precision)

        m_chk <- length(gene_idx_chk)
        ## The FUSED path also needs the pre-solve linear predictor for its rho
        ## accumulation, whose row sums run over the whole logical chunk, so it
        ## builds eta at chunk width; the default path builds it per sub-block.
        eta_chk <- if (fuse_rho && it > 1L) .eta_block(B, U, gene_idx_chk) else NULL
        if (fuse_rho) {
          ## num/den from the PRE-solve B,U (one-iteration lag vs the dense
          ## post-solve; same fixed point), with no convergence metric: the fused
          ## path uses a coefficient metric below.
          acc <- pace_rho_accumulate_cpp(
            if (it == 1L) empty_matrix else eta_chk, empty_matrix, Y, a_cache, gene_idx_chk[1L],
            offset_vec, add_rho, alpha[gene_idx_chk], mask_matrix, mask_index, num, den, disp_nb2,
            seed_iteration = (it == 1L), seed_previous = FALSE, n_cells = n, n_genes_in = m_chk,
            want_tail_counts = FALSE, n_threads = .pace_thread_count(n_threads))
          num <- acc$num
          den <- acc$den
          rm(acc)
        }
        ## The core walks the sub-blocks itself, writing each one straight into
        ## its own columns of a single z and w. The R loop this replaces allocated
        ## z_chk and w_chk per chunk, then allocated another pair per sub-block
        ## inside pace_working_response_cpp() and copied them in -- 1.25 GB of R
        ## allocation per chunk at 1.2M cells and a chunk of 64, before the
        ## per-sub-block temporaries. Sub-blocking still bounds the memory: eta is
        ## built sub_genes genes at a time unless the fused path already supplied
        ## it at chunk width, in which case a sub-block is a pointer into it.
        working <- pace_working_response_chunk_cpp(
          if (it == 1L || is.null(eta_chk)) empty_matrix else eta_chk,
          x1_or_empty, isTRUE(x1_is_unit), x_fixed_dense, p, B, Z, U,
          as.integer(gene_idx_chk), Y, a_cache, gene_idx_chk[1L], offset_vec, add_rho,
          alpha[gene_idx_chk], sample_weight_vec, disp_nb2, it == 1L, n,
          as.integer(sub_genes), .pace_thread_count(n_threads))
        z_chk    <- working$z
        w_chk    <- working$w
        colsum_w <- working$colsum_w
        rm(working, eta_chk)
        lam_chk <- lam_diag_mat[, gene_idx_chk, drop = FALSE]
        ## Keep the ridge out of single precision: below eps_float * sum(w) it is
        ## quantised away, at a threshold that falls as 1/n. The decision is made
        ## over the whole logical chunk, exactly as before the sub-block split.
        if (iter_precision != 0L &&
            .ridge_needs_double_from(suppressWarnings(min(lam_chk, na.rm = TRUE)),
                                     suppressWarnings(max(colsum_w))))
          iter_precision <- 0L

        per_gene_chk <- .solve_genes_chunk_multiblock(
          X_fixed, re$X_terms_list, re$cell_grp_list, re$cells_by_grp_list,
          w_chk, z_chk, lam_chk, re$blocks,
          gene_idx = seq_along(gene_idx_chk),
          n_threads = n_threads,
          interior_precision = iter_precision,
          BPPARAM = BPPARAM)
        ## LOSSLESS NaN GUARD (root-cause fix): the float (interior_precision=1)
        ## per-gene Cholesky in the compiled solve returns a NaN BLUP column when
        ## a gene's working-weight system is borderline in single precision (rare,
        ## run-to-run-variable because float OMP dgemm accumulation order is not
        ## bit-reproducible). A single NaN BLUP column makes Z%*%U[,g] NaN, which
        ## NaN-poisons den = rowSums(WA * a) for EVERY cell with nonzero ambient in
        ## that gene (i.e. nearly all cells) and would silently collapse add_rho to
        ## its prior panel-wide -- masquerading as "a few degenerate cells". We
        ## detect any non-finite-BLUP gene in this chunk and RE-SOLVE it in DOUBLE
        ## precision (interior_precision=0), which is numerically robust (the 1/tau
        ## ridge makes every per-gene system positive definite). Double re-solve was
        ## verified NaN-free under stress (lam up to 1e8, alpha up to 50). This keeps
        ## ALL genes valid -- no gene loses its BLUP -- and is a no-op when the float
        ## solve already returned finite values (the canonical case).
        ## A gene is bad if any of its beta or u entries is non-finite. colSums
        ## carries a NaN or an Inf straight through, so one pass over the two
        ## matrices answers it for every gene at once.
        is_bad <- function(fit) which(!is.finite(colSums(fit$B)) | !is.finite(colSums(fit$U)))
        bad_jj <- is_bad(per_gene_chk)
        if (length(bad_jj) && iter_precision != 0L) {
          if (verbose)
            cat(sprintf("    [nan-guard] it=%d chunk@%d: %d gene(s) NaN in float solve -> re-solving in double\n",
                        it, cs, length(bad_jj)))
          redo <- .solve_genes_chunk_multiblock(
            X_fixed, re$X_terms_list, re$cell_grp_list, re$cells_by_grp_list,
            w_chk[, bad_jj, drop = FALSE], z_chk[, bad_jj, drop = FALSE],
            lam_chk[, bad_jj, drop = FALSE], re$blocks,
            gene_idx = seq_along(bad_jj),
            n_threads = n_threads,
            interior_precision = 0L,
            BPPARAM = BPPARAM)
          per_gene_chk$B[, bad_jj] <- redo$B
          per_gene_chk$U[, bad_jj] <- redo$U
          per_gene_chk$Ainv_diag[, bad_jj] <- redo$Ainv_diag
          bad_jj <- is_bad(per_gene_chk)          ## did the repair actually take?
        }
        ## Anything still non-finite is unrecoverable and goes into the saved fit,
        ## so count it and say so on EVERY iteration -- including the last, where
        ## the interior is already double and the repair above cannot help. This
        ## silence is what let a fit report converged = TRUE with all-NaN genes.
        if (length(bad_jj)) {
          nan_genes_iter <- c(nan_genes_iter, gene_idx_chk[bad_jj])
          if (verbose)
            cat(sprintf("    [nan-guard] it=%d chunk@%d: %d gene(s) STILL non-finite after double solve\n",
                        it, cs, length(bad_jj)))
        }
        ## Whole blocks, not one gene at a time. The columns are already in gene
        ## order and contiguous, so these are three (or five) assignments instead
        ## of six per gene.
        random_rows <- (p + 1):(p + q)
        B[, gene_idx_chk]      <- per_gene_chk$B
        U[, gene_idx_chk]      <- per_gene_chk$U
        re_var[, gene_idx_chk] <- pmax(per_gene_chk$Ainv_diag[random_rows, , drop = FALSE], 0)
        if (last_iter) {
          se_B[, gene_idx_chk] <- sqrt(pmax(per_gene_chk$Ainv_diag[seq_len(p), , drop = FALSE], 0))
          se_U[, gene_idx_chk] <- sqrt(pmax(per_gene_chk$Ainv_diag[random_rows, , drop = FALSE], 0))
        }
        rm(per_gene_chk, z_chk, w_chk, lam_chk, colsum_w)
      }
    }

    ## ----- Per-cell contamination update (streaming rho accumulation) -----
    ## Dense computes (lines ~333-379), over the FULL n x G matrices:
    ##   mu_bio = exp(eta + offset); mu_spill = rho_i * a; mu = mu_bio + mu_spill;
    ##   wcnt   = NB precision at mu_tot = mu_bio + mu_spill (PREVIOUS rho);
    ##   WA     = wcnt * a  (masked to anchor genes);
    ##   num    = rowSums(WA * Y); den = rowSums(WA * a);
    ##   then EB-shrink rho. Because num / den are pure GENE-SUMS, accumulating
    ##   them chunk-by-chunk is exact (up to float summation order).
    ##
    ## NOTE: dense uses mu_spill computed from the rho of the PREVIOUS iteration
    ## (mu_spill is updated at the END of the loop, after this block uses it).
    ## We mirror that with mu_spill_chk = add_rho * a_chk using the current
    ## add_rho (which is still the previous iteration's value at this point).
    ## SPEED (pass-fusion): fold the rel_delta convergence diagnostic INTO this
    ## rho-accumulation chunk loop. eta_chk built here for current B,U is reused
    ## as the rel_delta "new" side (B/U are final for this iteration -- the alpha
    ## and tau updates do not change them -- so this is identical to recomputing
    ## eta from current B,U later). Only eta_prev (the denominator side) is
    ## recomputed: for iter 1 from the dense seed log(pmax(Y,0.5))-offset; for
    ## iter>1 from prev_B/prev_U. Metric is EXACTLY
    ##   max(|eta_new - eta_prev| / max(|eta_prev|, 1e-3)).
    ## num/den + convergence metric come from EITHER the fused Pass 1 above
    ## (fuse_rho=TRUE) OR this dedicated post-solve second pass (default). The
    ## accumulators were initialised before Pass 1.
    if (!fuse_rho) {
    ## One core call per logical chunk (pace::rho_accumulate): it reads the
    ## counts and the ambient field sparsely, rebuilds mu_tot from the PREVIOUS
    ## rho as the dense solver did, accumulates num/den in long double across
    ## the chunk's genes in column order, and returns the convergence metric
    ## max(|eta - prev_eta| / max(|prev_eta|, 1e-3)) over the same pass. A
    ## non-finite (cell, gene) contributes nothing to num/den, as rowSums with
    ## na.rm did, and is counted so a lossy fit cannot report convergence.
    ## The whole sweep runs in the core. Both eta matrices -- the current one and
    ## the previous iteration's -- are buffers reused across chunks rather than a
    ## pair of n x chunk R matrices built and dropped per chunk, and the running
    ## accumulators stay there instead of crossing back on every chunk.
    acc <- pace_rho_pass_cpp(
      x1_or_empty, isTRUE(x1_is_unit), x_fixed_dense, p, Z, B, U,
      if (prev_eta_set) prev_B else B, if (prev_eta_set) prev_U else U, prev_eta_set,
      Y, a_cache, offset_vec, add_rho, alpha, mask_matrix, mask_index, disp_nb2,
      RD_DIAG, n, as.integer(chunk_size), .pace_thread_count(n_threads))
    num <- acc$num
    den <- acc$den
    rel_delta <- max(rel_delta, acc$rel_delta_max)
    rd_n   <- rd_n   + acc$n_finite
    rd_sum <- rd_sum + acc$rel_delta_sum
    rd_nonfinite <- rd_nonfinite + acc$n_nonfinite
    if (RD_DIAG) {         ## tail-fraction counts: diagnostic only
      rd_g01 <- rd_g01 + acc$tail_counts[1L]
      rd_g05 <- rd_g05 + acc$tail_counts[2L]
      rd_g1  <- rd_g1  + acc$tail_counts[3L]
      rd_g10 <- rd_g10 + acc$tail_counts[4L]
    }
    rm(acc)
    } else {
      ## Fused path: num/den already accumulated in Pass 1. COEFFICIENT-based
      ## convergence metric -- the eta-based one is unusable here because the
      ## pre-solve eta equals the previous iteration's eta. rd = |U_new - U_prev|
      ## / max(|U_prev|, 1e-3) over the random-effect coefficients (they carry the
      ## cell x neighbour spatial signal; B is intercept-only).
      if (is.null(prev_U)) {                 ## iter 1: no previous coefficients yet
        rel_delta <- 1e3; rd_n <- length(U); rd_sum <- 1e3 * length(U)
      } else {
        dU   <- abs(U - prev_U) / pmax(abs(prev_U), 1e-3)
        finU <- is.finite(dU)
        rel_delta <- max(dU[finU]); rd_n <- sum(finU); rd_sum <- sum(dU[finU])
        if (RD_DIAG) {
          rd_g01 <- sum(dU[finU] > 0.01); rd_g05 <- sum(dU[finU] > 0.05)
          rd_g1  <- sum(dU[finU] > 0.1);  rd_g10 <- sum(dU[finU] > 1.0)
        }
        rm(dU, finU)
      }
    }
    ## rel_delta is now fully accumulated (Pass 2 or fused). After iter 1 the dense
    ## seed is no longer the previous eta (prev_B/prev_U are).
    prev_eta_set <- TRUE
    ## Empirical-Bayes shrink of the per-cell loading (pace::rho_shrink). A cell
    ## with non-finite accumulators carries no contamination information, so it
    ## gets the prior and drops out of the precision weights; the count is
    ## logged so a real divergence is visible rather than silently absorbed.
    shrunk  <- pace_rho_shrink_cpp(num, den)
    add_rho <- shrunk$rho
    n_bad   <- shrunk$n_nonfinite
    if (verbose && n_bad > 0L)
      cat(sprintf("    [percell_bleed] it=%d guarded %d/%d non-finite den/num cells (rho=prior there)\n",
                  it, n_bad, length(den)))
    if (verbose && it <= 3) {
      ## Streamed contam_frac diagnostic (matches dense print).
      spill_sum <- numeric(n); tot_sum <- numeric(n)
      for (cs in chk_starts) {
        gene_idx_chk <- cs:min(cs + chunk_size - 1L, g_n)
        eta_chk <- .eta_block(B, U, gene_idx_chk)
        pass <- pace_final_pass_statistics_cpp(eta_chk, offset_vec, a_cache, gene_idx_chk[1L],
                                               add_rho, diagnostic_group, 1L,
                                               want_groups = FALSE, return_matrices = FALSE)
        spill_sum <- spill_sum + pass$spill_row_sum
        tot_sum   <- tot_sum   + pass$total_row_sum
        rm(eta_chk, pass)
      }
      fr <- spill_sum / pmax(tot_sum, 1e-9)
      cat(sprintf("    [percell_bleed] it=%d  rho_i [%.4f,%.4f] med=%.4f  contam_frac med=%.3f q90=%.3f\n",
                  it, min(add_rho), max(add_rho), stats::median(add_rho),
                  stats::median(fr), stats::quantile(fr, 0.9)))
    }
    rm(num, den)

    ## ----- alpha MLE per gene (streamed mu reconstruction per gene) -----
    ## Dense computes alpha on the FULL mu = mu_bio + rho_i*a (the NEW add_rho).
    ## We reconstruct mu[, gi] per chunk and feed each gene column to the MLE.
    ## SPEED 4 (guarded): the alpha MLE is ~37% of per-iter cost and converges
    ## quickly, so we only re-fit it for the first `alpha_warmup` iterations (and
    ## the last iteration); otherwise alpha / prev_alpha are left unchanged.
    ## alpha_warmup = Inf disables the freeze (always update = exact).
    update_alpha <- (it <= alpha_warmup) || last_iter
    if (update_alpha) {
    ## The whole sweep runs in the core: each gene's fitted mean is rebuilt from
    ## its own column and its dispersion minimised over log alpha, on the core's
    ## threads, and eta is a buffer reused across chunks rather than an n x chunk
    ## R matrix built and dropped per chunk.
    fitted <- pace_dispersion_pass_cpp(
      x1_or_empty, isTRUE(x1_is_unit), x_fixed_dense, p, Z, B, U,
      Y, a_cache, offset_vec, add_rho, disp_nb2, alpha_zero_collapse, alpha_max_n,
      alpha_fast_density, n, as.integer(chunk_size), .pace_thread_count(n_threads))
    if (fitted$n_noninteger > 0)
      warning(sprintf("iter %d: %.0f gene(s) have counts that are not whole numbers; their ",
                      it, fitted$n_noninteger),
              "dispersion is undefined and keeps its previous value.", call. = FALSE)
    alpha_new <- fitted$alpha
    rm(fitted)
    alpha <- alpha_new
    alpha[!is.finite(alpha)] <- prev_alpha[!is.finite(alpha)]
    alpha <- pmin(pmax(alpha, 1e-4), 50)
    prev_alpha <- alpha
    } else if (verbose) {
      cat(sprintf("    [alpha] it=%d > alpha_warmup=%s: alpha FROZEN\n",
                  it, format(alpha_warmup)))
    }

    ## ----- tau update (identical to dense) -----
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
                           d0_min = d0_min_env, p_fixed = p)
        } else {
          .adaptive_tau_half_cauchy(s2_b, blk_i$K_terms, blk_i$K_groups)
        }
        tau_g_array[rng, ] <- tau_b_g
      }
    }
    if (!is.null(data_informed_W)) {
      stopifnot(identical(dim(data_informed_W), dim(tau_g_array)))
      tau_g_array <- tau_g_array * data_informed_W
      tau_g_array <- pmax(tau_g_array, 1e-8)
    }
    ## Every tau write is a pmax against a floor; without a ceiling a column
    ## identified only by the ridge climbs by mean(u^2) per iteration, and the
    ## inflated posterior variance feeds the next update.
    tau_g_array <- .clamp_tau(tau_g_array, tau_max, it)
    hist$tau_max_seen[it] <- max(tau_g_array, na.rm = TRUE)

    ## ----- rel_delta convergence diagnostic -----
    ## rel_delta was computed (fused) inside the rho-accumulation loop above,
    ## reusing the current-B,U eta_chk built there. Same exact metric:
    ##   max(|eta_new - prev_eta| / max(|prev_eta|, 1e-3)). Stash prev_B/prev_U
    ## for the NEXT iteration's denominator.
    prev_B <- B; prev_U <- U

    hist$tau_blocks[[it]] <- tau_blocks
    hist$alpha[[it]]      <- alpha
    hist$rel_delta[it]    <- rel_delta                               ## L-inf max (diagnostic)
    hist$rel_delta_mean[it] <- if (rd_n > 0) rd_sum / rd_n else rel_delta  ## mean (early-stop metric)
    nan_genes_iter <- unique(nan_genes_iter)
    hist$n_nan_genes[it]           <- length(nan_genes_iter)
    hist$n_nonfinite[it]           <- rd_nonfinite
    if (length(nan_genes_iter)) {
      warning(sprintf(
        "iter %d: %d gene(s) carry non-finite coefficients (%s%s); the fit will not be reported as converged.",
        it, length(nan_genes_iter),
        paste(utils::head(colnames(Y)[nan_genes_iter], 5), collapse = ", "),
        if (length(nan_genes_iter) > 5) ", ..." else ""), call. = FALSE)
    }
    if (RD_DIAG && rd_n > 0)
      cat(sprintf("  [rd-diag] it=%d  mean=%.4g  max=%.3g  frac>0.01=%.2e  >0.05=%.2e  >0.1=%.2e  >1=%.2e  (N=%.2e)\n",
                  it, rd_sum / rd_n, rel_delta, rd_g01/rd_n, rd_g05/rd_n, rd_g1/rd_n, rd_g10/rd_n, rd_n))
    if (verbose) {
      ## Report the max as well as the median: a runaway confined to the
      ## intercept row does not move a median taken over the whole block.
      tau_med <- vapply(tau_blocks, function(m) stats::median(m), numeric(1))
      cat(sprintf("  [mvpql.streaming] iter %d  rel_delta[mean]=%.3g (max=%.3g)  alpha[med]=%.2f  tau=[%s] (max %.3g)  (%.1fs)\n",
                  it, hist$rel_delta_mean[it], rel_delta, median(alpha),
                  paste(sprintf("%.3f", tau_med), collapse = ","),
                  hist$tau_max_seen[it],
                  as.numeric(difftime(Sys.time(), t_it, units = "secs"))))
    }
    ## SPEED 4 (guarded): break once this iteration was flagged as the early-stop
    ## last_iter (SEs + final alpha already computed above this iteration).
    if (stop_now) {
      converged <- TRUE
      if (verbose)
        cat(sprintf("  [mvpql.streaming] EARLY STOP at iter %d (mean rel_delta[%d]=%.3g < %.3g; L-inf max=%.3g)\n",
                    it, it - 1L, hist$rel_delta_mean[it - 1L], early_stop_tol, hist$rel_delta[it - 1L]))
      break
    }
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
