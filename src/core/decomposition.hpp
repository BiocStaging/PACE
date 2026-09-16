// decomposition.hpp -- the arithmetic of the PACE reporting layer: the per-gene
// variance decomposition, its aggregate tables, the observed single-frame shares,
// the per-pair driver scores and the Pratt pair-variance attribution.
//
// Every function here works on per-cell-type statistics and on the fit's own
// parameter matrices, never on a cells x genes matrix.
#ifndef PACE_DECOMPOSITION_HPP
#define PACE_DECOMPOSITION_HPP

#include <cstdint>

#include "core_types.hpp"
#include "count_stats.hpp"

namespace pace {

// Inputs of the per-(focal, gene) variance decomposition.
//
// For focal cell type c (its n_c cells) and gene g, with the kernel columns
// N_i (a K-vector per cell), the slope BLUPs u_cg and their standard errors:
//
//   V_state_baseline  = var_i(N_i' u_cg) + sum_k se_kcg^2 var_i(N_ik)
//                     = u_cg' S_c u_cg + diag(S_c)' se_cg^2,  S_c = Cov_{i in c}(N_i)
//   V_state_responder = r_cg' S^R_c r_cg + diag(S^R_c)' se_rcg^2
//   V_spill           = toff_var[c, g]                      (per-cell contamination)
//                     = b_g' Cov_{i in c}(X_spill,i) b_g     (legacy _near fixed effects)
//   V_disp            = log(1 + (1 + max(alpha_g, 0)) / max(mu_mean[c, g], 1e-9))   (NB1)
//                     = log(1 + 1 / max(mu_mean[c, g], 1e-9) + max(alpha_g, 0))     (NB2)
//   celltype_offset_sq = u_intercept,cg^2 + se_intercept,cg^2
//
// (var(a' x) = a' Cov(x) a for a fixed vector a, which is why no per-cell value
// is needed.) The five percentages divide each block by their sum.
//
// Index conventions: `focal_group`, `slope_rows`, `responder_rows`,
// `responder_keep` and `intercept_rows` are 1-based, as R indexes; a
// non-positive entry means "absent". Matrices are column-major.
struct DecompositionInput {
  std::int64_t n_genes = 0;
  int n_groups = 0;                  // cell types with count means and covariance slices
  int n_focals = 0;                  // reported focal cell types
  Span<const double> ct_means;       // n_groups x n_genes: mean count per type and gene
  Span<const int> group_size;        // n_groups: cells per type (0: means read as 0)
  Span<const int> focal_group;       // n_focals: the group each focal refers to
  Span<const int> n_focal;           // n_focals: cells of each focal
  Span<const double> u;              // q x n_genes: the BLUPs (fit$U)
  Span<const double> se_u;           // q x n_genes: their standard errors
  std::int64_t q = 0;
  std::int64_t n_kernel = 0;
  Span<const int> slope_rows;        // n_kernel x n_focals: rows of u holding the slopes
  Span<const double> kernel_cov;     // n_groups slices of n_kernel x n_kernel
  std::int64_t n_responder = 0;      // 0 when the fit has no condition interaction
  Span<const int> responder_rows;    // n_responder x n_focals
  Span<const int> responder_keep;    // n_responder: which kernel columns they match
  Span<const double> responder_cov;  // n_groups slices of n_kernel x n_kernel
  Span<const int> intercept_rows;    // n_focals: row of u holding the type offset
  Span<const double> toff_var;       // n_focals x n_genes, empty unless used
  std::int64_t n_spill = 0;          // legacy spillover fixed effects
  Span<const double> spill_cov;      // n_groups slices of n_spill x n_spill
  Span<const double> beta_spill;     // n_spill x n_genes
  Span<const double> mu_mean;        // n_focals x n_genes: mean fitted mean
  Span<const double> alpha;          // n_genes: the NB dispersion
  bool nb1 = true;
};

// Outputs, each n_focals * n_genes with the row of focal f and gene j at
// f * n_genes + j, the order the reported table uses. `keep` is 0 for a focal
// with fewer than five cells and for a gene whose blocks do not sum to a finite
// positive total, exactly as the R implementation dropped those rows.
struct DecompositionOutput {
  Span<double> v_state_baseline;
  Span<double> v_state_responder;
  Span<double> v_spill;
  Span<double> v_disp;
  Span<double> celltype_offset_sq;
  Span<double> focal_mean;
  Span<double> max_other_mean;
  Span<double> spec;
  Span<double> focal_other_ratio;
  Span<int> is_contaminated;
  Span<double> pct_celltype;
  Span<double> pct_state;
  Span<double> pct_responder;
  Span<double> pct_spill;
  Span<double> pct_residual;
  Span<int> keep;
};

Status variance_decomposition(const DecompositionInput& input, const DecompositionOutput& output,
                              int n_threads, const InterruptCheck& interrupted);

// The aggregate tables of the decomposition, over the rows of one focal:
//   mean    = mean(block %, na.rm = TRUE)                       per block
//   specw   = sum(block % * spec^2) / max(sum(spec^2), 1e-12)   per block
//   pooled  = 100 * sum(block) / max(sum(Total), 1e-12)         per block
//   Total   = V_state_baseline + V_state_responder + V_spill + V_disp + celltype_offset_sq
// with every sum in long double and NA dropped, as R's sum(na.rm = TRUE).
//
// Shapes: `focal_code` has n_rows entries in [0, n_focals); `pct` and
// `components` are n_rows x 5 column-major, in the block order
// (cell type, spatial state, responder spatial state, spillover, residual) for
// `pct` and (celltype_offset_sq, V_state_baseline, V_state_responder, V_spill,
// V_disp) for `components`; the three tables are n_focals x 5; `total` is n_rows;
// `total_ss`, `n_rows_out` and `n_specific` are n_focals. `specw` may be empty.
Status decomposition_aggregates(Span<const int> focal_code, int n_focals, std::int64_t n_rows,
                                Span<const double> pct, Span<const double> components,
                                Span<const double> spec, Span<const int> is_contaminated,
                                Span<double> mean_table, Span<double> specw_table,
                                Span<double> pooled_table, Span<double> total,
                                Span<double> total_ss, Span<int> n_rows_out,
                                Span<int> n_specific);

// The four-block view: the responder spatial state folded into the spatial block.
//   Total   = celltype_offset_sq + V_state_baseline + V_state_responder + V_spill + V_disp
//   block % = 100 * block / max(Total, 1e-12)
// Shapes: every span has n_rows entries.
Status four_block_shares(Span<const double> celltype_offset_sq, Span<const double> v_state_baseline,
                         Span<const double> v_state_responder, Span<const double> v_spill,
                         Span<const double> v_disp, std::int64_t n_rows, Span<double> total,
                         Span<double> v_state, Span<double> pct_celltype, Span<double> pct_state,
                         Span<double> pct_spill, Span<double> pct_residual);

// The observed single-frame decomposition, for focal type c and gene g, with
// y = log1p(CP10k) and the statistics of single_frame_statistics():
//
//   SS_lineage = n_c (focal_mean - global_mean)^2
//   denom      = SS_lineage + SS_within
//   Cell type %= 100 SS_lineage / denom
//   block %    = 100 SS_within p_block / denom,  p_block = V_block / sum of the fit's V
//
// so the magnitudes are observed and only the within-cell-type split comes from
// the fit. A type with fewer than five cells gets keep = 0.
//
// Shapes: focal_mean, within_ss and the four V spans are n_groups x n_genes
// column-major (the V spans aligned to the same genes, NaN where the fit has no
// row); global_mean and the outputs follow the same n_groups x n_genes layout,
// except `keep`, which has n_groups entries. `v_responder` may be empty.
Status single_frame_shares(Span<const double> focal_mean, Span<const double> within_ss,
                           Span<const double> global_mean, Span<const int> group_size,
                           int n_groups, std::int64_t n_genes, Span<const double> v_state,
                           Span<const double> v_responder, Span<const double> v_spill,
                           Span<const double> v_disp, Span<double> pct_celltype,
                           Span<double> pct_spatial, Span<double> pct_responder,
                           Span<double> pct_spill, Span<double> pct_residual,
                           Span<double> ss_lineage, Span<double> ss_within, Span<double> denom,
                           Span<int> keep, int n_threads, const InterruptCheck& interrupted);

// SS-weighted per-focal blocks of the observed single frame (Goldstein pooling):
//   ss_block = block % / 100 * denom            per row
//   block    = 100 * sum(ss_block) / sum(denom) per focal
// The sums keep NA, as the R implementation's sum() did.
// Shapes: `focal_code` has n_rows entries in [0, n_focals); `pct` is
// n_rows x n_blocks column-major; `denom` has n_rows entries; `table` is
// n_focals x n_blocks.
Status single_frame_focal_blocks(Span<const int> focal_code, int n_focals, std::int64_t n_rows,
                                 Span<const double> pct, int n_blocks, Span<const double> denom,
                                 Span<double> table);

// Per-pair driver scores for one (focal, neighbour) pair, over the pair's genes:
//   MCSD    = b^2 spec^2 max(focal_mean, 0)         (the locked driver score)
//   MCSD4   = b^2 spec^4 max(focal_mean, 0)
//   V_resid = log(1 + (1 + max(alpha, 0)) / max(mu_bar, 1e-6))
//   V_S     = b^2 var_n                              (no condition), or u_raw^2 var_n
//   V_RxS   = b^2 var_rn                             (condition fits only)
//   R2      = V / max(V_total, 1e-12)
// with b the shrunken slope and u_raw the unshrunken BLUP. `u_raw`, `v_rxs` and
// `r2_rxs` are used only when has_responder; every span has n entries.
Status driver_scores(Span<const double> estimate_shrunk, Span<const double> spec,
                     Span<const double> focal_mean, Span<const double> mu_bar,
                     Span<const double> alpha, Span<const double> u_raw, double var_n,
                     double var_rn, bool has_responder, std::int64_t n, Span<double> mcsd,
                     Span<double> mcsd4, Span<double> v_resid, Span<double> v_s,
                     Span<double> v_rxs, Span<double> v_total, Span<double> r2_s,
                     Span<double> r2_rxs);

// Sample covariance of selected columns of a sparse cells x terms matrix over
// selected rows, with stats::cov semantics (see group_covariances). `scale`, when
// not empty, multiplies every selected column elementwise by one value per
// selected row, so var(R N) is one call. Fewer than two rows gives NaN.
//
// Shapes: `rows` and `cols` are 1-based indices into `matrix`; `rows` must be
// strictly increasing, so each column is gathered by one walk of its stored
// entries and nothing of the matrix's height is allocated. `scale` is empty or
// one value per selected row; `covariance` is n_cols * n_cols column-major.
Status subset_covariance(const CscView& matrix, Span<const int> rows, Span<const int> cols,
                         Span<const double> scale, Span<double> covariance);

// The values of one column of a sparse matrix at selected, strictly increasing
// (1-based) rows.
Status subset_column(const CscView& matrix, Span<const int> rows, int column, Span<double> values);

// The 1-based rows at which one column of a sparse matrix is not zero.
Status column_nonzero_rows(const CscView& matrix, int column, Span<int> rows, std::int64_t* n_found);

// Pratt's signed pair-level attribution of one focal's spatial variance
// (Pratt 1987): with Sigma the covariance of the focal's kernel columns and U the
// slope BLUPs (non-finite entries read as 0),
//   V_pair[k]   = sum_g U[k, g] (Sigma U)[k, g]        (signed, sums to the block)
//   V_diag[k]   = sum_g U[k, g]^2 Sigma[k, k]          (the uncorrelated part)
//   share %     = 100 V_pair[k] / sum_k V_pair[k]
//   cross_cov % = 100 (sum V_pair - sum V_diag) / sum V_pair
// Shapes: sigma is n_pairs x n_pairs, u is n_pairs x n_genes, the outputs have
// n_pairs entries.
Status pair_variance_pratt(Span<const double> sigma, Span<const double> u, std::int64_t n_pairs,
                           std::int64_t n_genes, Span<double> v_pair, Span<double> v_diag,
                           Span<double> share_pct, Span<double> diag_pct, Span<double> sign,
                           double* v_total, double* v_diag_total, double* cross_cov_pct,
                           int* n_negative);

}  // namespace pace

#endif  // PACE_DECOMPOSITION_HPP
