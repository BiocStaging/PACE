// statistics.hpp -- the per-cell-type statistics the reporting layers read from
// a fit, accumulated without holding an n x G matrix in R.
#ifndef PACE_STATISTICS_HPP
#define PACE_STATISTICS_HPP

#include <cstdint>

#include "core_types.hpp"

namespace pace {

// True when every value is finite (no NA, NaN or Inf). Used to validate large
// vectors without materialising R's logical copy of them.
bool all_finite(Span<const double> values);

// Per-group mean and sample variance of every column of a dense n x p matrix,
// exactly as R computes them on the group's rows:
//
//   mean_gj = mean(x[group == g, j], na.rm = TRUE)
//           = s / m refined once by s += sum_k (x_kj - s), both in long double
//   var_gj  = stats::var(x[group == g, j], na.rm = TRUE)
//           = sum_k (x_kj - mean_gj)^2 / (m - 1), the square in double and the
//             sum in long double, over the group's rows in ascending order
//
// with m the number of values that are not NA or NaN (na.rm drops them, exactly
// as R does, before either statistic is computed). `n_used` reports m per slot,
// so the caller can render R's NA for m < 2 (var of fewer than two values).
//
// An empty `values` span means an all-zero matrix, so a caller with no technical
// offset need not allocate one.
//
// Shapes: values is n * p column-major or empty; group has n codes in
// [-1, n_groups) (-1: no group); `mean`, `variance` and `n_used` are each either
// empty (not wanted) or n_groups * p, column-major (groups x columns).
// `any_nonzero`, when not null, reports whether any entry of the whole matrix is
// non-zero, ignoring NA and NaN, as R's any(x != 0, na.rm = TRUE).
Status dense_group_moments(Span<const double> values, std::int64_t n, std::int64_t p,
                           Span<const int> group, int n_groups, Span<double> mean,
                           Span<double> variance, Span<int> n_used, bool* any_nonzero,
                           int n_threads, const InterruptCheck& interrupted);

// One gene chunk of the streaming solver's final pass: build the fitted means and
// the contamination log-offset and accumulate everything the fit keeps, without
// materialising the four n x n_genes working matrices R used to hold.
//
// For cell i and gene j of the chunk (eta = X B + Z U, supplied by the caller):
//   mu_bio  = max(exp(eta_ij + offset_i), 1e-6)
//   spill   = max(ambient_ij * rho_i, 0)
//   mu      = max(mu_bio + spill, 1e-6)
//   toff    = log1p(spill / mu_bio)
// and then, matching the R accumulation term for term:
//   mu_group_sum[g, j]  = colSums(mu[group == g, ])     (long double, ascending rows)
//   toff_variance[g, j] = stats::var(toff[group == g, j])
//   mu_column_sum[j]    = colSums(mu)
//   spill_row_sum[i]    = rowSums(spill)   over this chunk's genes
//   total_row_sum[i]    = rowSums(mu)      over this chunk's genes
//   any_nonzero         = any(toff != 0)
//
// Deliberately serial: the row sums accumulate across the chunk's genes in
// column order, so a thread split would change their last bits with the thread
// count. The work is one pass over n x n_genes, the same pass R made.
//
// Shapes: eta and ambient are n * n_genes column-major; offset, rho and group
// have n entries; mu_group_sum and toff_variance are n_groups * n_genes;
// mu_column_sum has n_genes; spill_row_sum and total_row_sum have n; mu_out and
// toff_out are either empty (not wanted) or n * n_genes. A group with fewer than
// two cells gets a NaN variance, which the binding reports as R's NA.
// Non-finite inputs are refused.
Status final_pass_statistics(Span<const double> eta, std::int64_t n, std::int64_t n_genes,
                             Span<const double> offset, Span<const double> ambient,
                             Span<const double> rho, Span<const int> group, int n_groups,
                             Span<double> mu_group_sum, Span<double> toff_variance,
                             Span<double> mu_column_sum, Span<double> spill_row_sum,
                             Span<double> total_row_sum, Span<double> mu_out,
                             Span<double> toff_out, bool* any_nonzero,
                             const InterruptCheck& interrupted);

}  // namespace pace

#endif  // PACE_STATISTICS_HPP
