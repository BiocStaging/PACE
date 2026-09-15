// count_stats.hpp -- per-group column statistics from sparse counts.
#ifndef PACE_COUNT_STATS_HPP
#define PACE_COUNT_STATS_HPP

#include <cstdint>

#include "core_types.hpp"

namespace pace {

// Compressed-column view of a cells x genes count matrix (the dgCMatrix slots):
// column_pointer has n_genes + 1 entries; row_index and values have nnz entries;
// row indices are 0-based and ascending within each column.
struct CscView {
  Span<const int> column_pointer;
  Span<const int> row_index;
  Span<const double> values;
  std::int64_t n_rows = 0;
  std::int64_t n_cols = 0;
};

// Mean of every column within each cell group, as R's
//   colMeans(Y[group == g, , drop = FALSE])          (detection = false)
//   colMeans(Y[group == g, , drop = FALSE] > 0)      (detection = true)
// The column sum over the group's cells is accumulated in long double in row
// order (zeros add nothing), then divided by the group size in long double and
// rounded to double, which is exactly what colMeans() does. A group with no
// cells gives NaN (0/0), as colMeans() of a zero-row matrix. Cells with group -1
// belong to no group. Non-finite values (NA, NaN, Inf) are refused.
//
// Shapes: group has n_rows entries with codes in [-1, n_groups). `means` is
// caller-owned with n_groups * n_cols entries, column-major (groups x genes).
Status group_column_means(const CscView& counts, Span<const int> group, int n_groups,
                          bool detection, Span<double> means, int n_threads,
                          const InterruptCheck& interrupted);

// Per-group sample covariance of the columns of a dense n x p matrix, as R's
// stats::cov(X[group == g, , drop = FALSE]) computes it (src/library/stats/src/cov.c):
//
//   mean_j = s / n_g with s the long double column sum, refined once by
//            mean_j += sum_k (x_kj - mean_j) / n_g (long double);
//   cov_jl = [ sum_k (x_kj - mean_j) (x_kl - mean_l) ] / (n_g - 1),
//
// the product in double and the sum in long double, over the group's rows in
// order. A group with fewer than 2 rows gives NaN. Non-finite values are refused.
//
// Shapes: values is n * p column-major (caller-owned); group has n codes in
// [-1, n_groups); covariance is caller-owned with n_groups * p * p entries, the
// p x p matrix of group g stored column-major at offset g * p * p.
Status group_covariances(Span<const double> values, std::int64_t n, std::int64_t p,
                         Span<const int> group, int n_groups, Span<double> covariance,
                         int n_threads, const InterruptCheck& interrupted);

// Observed single-frame statistics of a cells x genes count matrix:
//
//   y_ig = log1p(count_ig * (1e4 / n_count_i))          (log1p CP10k; 0 where count is 0)
//   focal_mean[g, j]  = (long double sum_{i in g} y_ij) / n_g
//   within_ss[g, j]   = sum_{i in g} (y_ij - focal_mean[g, j])^2
//                     = sum over the group's stored entries of (y - m)^2 + (n_g - stored) * m^2
//   global_mean[j]    = (long double sum_i y_ij) / n
//
// n_count must be finite and > 0 for every cell (a zero library size makes y
// undefined). Groups with no cells give NaN means and 0 sums of squares.
// Shapes: n_count and group have n entries; focal_mean and within_ss have
// n_groups * n_genes entries (column-major groups x genes); global_mean has n_genes.
Status single_frame_statistics(const CscView& counts, Span<const double> n_count,
                               Span<const int> group, int n_groups, Span<double> focal_mean,
                               Span<double> within_ss, Span<double> global_mean, int n_threads,
                               const InterruptCheck& interrupted);

}  // namespace pace

#endif  // PACE_COUNT_STATS_HPP
