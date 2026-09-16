// irls_chunk.hpp -- the per-chunk arithmetic of the streaming PQL solver.
//
// The solver walks the genes in logical chunks. Everything each chunk needs is
// built here from the counts and the ambient field read SPARSELY (one gene
// column at a time, expanded into a scratch vector), so the n x chunk working
// matrices R used to hold -- the counts, the ambient block, the two parts of the
// fitted mean, the contamination offset -- are never formed.
#ifndef PACE_IRLS_CHUNK_HPP
#define PACE_IRLS_CHUNK_HPP

#include <cstdint>

#include "core_types.hpp"
#include "count_stats.hpp"

namespace pace {

// The working response and weights of one gene block, as the dense solver
// formed them:
//
//   seed iteration:  mu_bio = max(y, 0.5),  mu = mu_bio
//   later:           mu_bio = max(exp(eta + offset), 1e-6)
//                    spill  = max(ambient * rho, 0)
//                    mu     = max(mu_bio + spill, 1e-6)
//   z = (log(mu_bio) - offset) + (y - mu) / mu_bio
//   w = (mu_bio^2 / mu) / (1 + alpha)          (NB1)
//     = mu_bio^2 / (mu (1 + mu alpha))         (NB2)
//   w = w * sample_weight                      (when given)
//
// The two w expressions keep the association the R code used, which is not the
// same for the two families. `colsum_w[j]` is sum_i w[i, j] in long double, the
// statistic the float/double gate reads; it is per gene, so it is unaffected by
// how the caller splits the chunk.
//
// Shapes: `eta` is n * n_genes column-major, or empty at the seed iteration,
// where `ambient` is ignored; `offset` and `rho` have n entries; `alpha` has
// n_genes; `sample_weight` has n entries or is empty; `z` and `w` are
// n * n_genes; `colsum_w` has n_genes. Non-finite inputs are refused.
Status working_response(Span<const double> eta, const GeneBlock& counts,
                        const GeneBlock& ambient, Span<const double> offset,
                        Span<const double> rho, Span<const double> alpha,
                        Span<const double> sample_weight, bool nb2, bool seed_iteration,
                        std::int64_t n, std::int64_t n_genes, Span<double> z, Span<double> w,
                        Span<double> colsum_w, int n_threads, const InterruptCheck& interrupted);

// The per-cell contamination accumulators and the convergence metric of one
// logical chunk, as the dense solver's second pass computed them:
//
//   mu_bio = max(exp(eta + offset), 1e-6)
//   spill  = max(ambient * rho, 0)                      (the PREVIOUS rho)
//   mu_tot = max(mu_bio + spill, 1e-8)                  (max(y, 0.5) at the seed)
//   wcnt   = (1 / mu_tot) / (1 + alpha)                 (NB1)
//          = 1 / (mu_tot (1 + mu_tot alpha))            (NB2)
//   WA     = wcnt * ambient, zeroed outside the anchor mask
//   num   += rowSums(WA * y),  den += rowSums(WA * ambient)     (NA dropped)
//
//   rd        = |eta - prev_eta| / max(|prev_eta|, 1e-3)
//   rel_delta = max(rd), and the mean over the finite entries drives the early stop
//
// The row sums accumulate in long double across the whole chunk in column
// order, which is what R's rowSums() did, so one call must cover one logical
// chunk. `num` and `den` are added to, not overwritten.
//
// The anchor mask is either per cell type (`mask` is n_types x n_genes with
// `mask_index` giving each cell's 1-based type) or per cell (n x n_genes with
// an empty `mask_index`); an empty `mask` means no masking.
//
// `prev_eta` is n * n_genes, or empty: with `seed_previous` the previous
// linear predictor is the solver's own seed, log(max(y, 0.5)) - offset, and
// without it the convergence metric is not computed at all (the fused path uses
// a coefficient metric instead).
//
// Shapes: `eta` is n * n_genes; `num` and `den` have n entries;
// `rel_delta_max`, `rel_delta_sum`, `n_finite` and `n_nonfinite` are single
// accumulators the caller owns; `tail_counts` is empty or four counters, the
// entries above 0.01, 0.05, 0.1 and 1 (the R_RD_DIAG diagnostic).
Status rho_accumulate(Span<const double> eta, Span<const double> prev_eta,
                      const GeneBlock& counts, const GeneBlock& ambient,
                      Span<const double> offset, Span<const double> rho,
                      Span<const double> alpha, Span<const double> mask,
                      Span<const int> mask_index, std::int64_t n_mask_rows, bool nb2,
                      bool seed_iteration, bool seed_previous, std::int64_t n,
                      std::int64_t n_genes, Span<double> num,
                      Span<double> den, double* rel_delta_max, double* rel_delta_sum,
                      std::int64_t* n_finite, std::int64_t* n_nonfinite, Span<double> tail_counts,
                      const InterruptCheck& interrupted);

// The fitted mean of one gene, mu = max(max(exp(eta + offset), 1e-6) +
// max(ambient rho, 0), 1e-6), and that gene's counts, for the dispersion MLE.
// `eta_column` and the outputs have n entries.
Status fitted_mean_column(Span<const double> eta_column, const GeneBlock& counts,
                          const GeneBlock& ambient, Span<const double> offset,
                          Span<const double> rho, std::int64_t n, std::int64_t gene,
                          Span<double> mu, Span<double> y);

// The empirical-Bayes shrink of the per-cell contamination loading, as the
// dense solver applied it:
//   ratio  = num / den where both are finite and den > 1e-12, else 0
//   raw    = max(ratio, 0)
//   weight = den where finite and den > 0, else 0
//   bar    = sum(weight raw) / max(sum(weight), 1e-12)
//   rho    = max((weight raw + den0 bar) / (weight + den0), 0)
// with den0 the 10th percentile of the positive weights (R's type-7 quantile).
// With no positive weight the loading is zero everywhere. `n_nonfinite` reports
// the cells whose accumulators were not finite, which the caller logs.
Status rho_shrink(Span<const double> num, Span<const double> den, std::int64_t n,
                  Span<double> rho, std::int64_t* n_nonfinite, double* rho_bar,
                  double* den_quantile);

}  // namespace pace

#endif  // PACE_IRLS_CHUNK_HPP
