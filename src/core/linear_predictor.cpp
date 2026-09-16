// fp_no_contract.hpp must come first: it forbids FMA contraction for the rest of
// the translation unit, which is what keeps this bit-identical to the R it
// replaces on FMA targets.
#include "fp_no_contract.hpp"

#include "linear_predictor.hpp"

#include <cstring>

#include "thread_pool.hpp"

namespace pace {

Status eta_block(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                 std::int64_t p, Span<const double> b, const CscView& z,
                 Span<const double> u, Span<const int> genes, std::int64_t n,
                 std::int64_t n_genes_total, Span<double> eta, int n_threads,
                 const InterruptCheck& interrupted) {
  const std::int64_t n_chunk = genes.size;
  if (n < 0 || p < 1)
    return Status::failure(StatusCode::invalid_argument, "eta_block: n and p must be positive.");
  if (z.n_rows != n)
    return Status::failure(StatusCode::invalid_argument,
                           "eta_block: Z has a different number of rows than n.");
  const std::int64_t q = z.n_cols;
  if (eta.size != n * n_chunk)
    return Status::failure(StatusCode::invalid_argument,
                           "eta_block: eta is not n * length(genes).");
  if (b.size != p * n_genes_total)
    return Status::failure(StatusCode::invalid_argument,
                           "eta_block: B is not p * n_genes_total.");
  if (u.size != q * n_genes_total)
    return Status::failure(StatusCode::invalid_argument,
                           "eta_block: U is not q * n_genes_total.");
  if (p == 1) {
    if (!x1_is_unit && x1.size != n)
      return Status::failure(StatusCode::invalid_argument,
                             "eta_block: x1 must have n entries when p == 1.");
  } else if (x_fixed.size != n * p) {
    return Status::failure(StatusCode::invalid_argument,
                           "eta_block: X_fixed is not n * p.");
  }
  for (std::int64_t j = 0; j < n_chunk; ++j) {
    const int gene = genes[j];
    if (gene < 0 || gene >= n_genes_total)
      return Status::failure(StatusCode::invalid_argument,
                             "eta_block: a gene index is out of range.");
  }
  if (n_chunk == 0) return Status::success();

  const Status status = parallel_for(
      n_chunk, n_threads, 1,
      [&](std::int64_t first, std::int64_t last) {
        for (std::int64_t j = first; j < last; ++j) {
          const std::int64_t gene = genes[j];
          double* column = eta.data + j * n;

          // ---- the Z part, accumulated in Z's own column order -------------
          // A zero coefficient contributes exactly zero to a finite sum, so the
          // column is skipped rather than walked; the core refuses non-finite
          // input elsewhere, so 0 * z can never be NaN here.
          std::memset(column, 0, static_cast<std::size_t>(n) * sizeof(double));
          const double* u_gene = u.data + gene * q;
          for (std::int64_t k = 0; k < q; ++k) {
            const double coefficient = u_gene[k];
            if (coefficient == 0.0) continue;
            const int begin = z.column_pointer[k];
            const int end = z.column_pointer[k + 1];
            for (int index = begin; index < end; ++index) {
              column[z.row_index[index]] += z.values[index] * coefficient;
            }
          }

          // ---- the fixed part, added once the Z sum is complete -------------
          if (p == 1) {
            const double coefficient = b.data[gene];
            if (x1_is_unit) {
              for (std::int64_t i = 0; i < n; ++i) column[i] += coefficient;
            } else {
              for (std::int64_t i = 0; i < n; ++i) column[i] += x1.data[i] * coefficient;
            }
          } else {
            // Each element's fixed part is summed over p on its own before it
            // reaches the accumulator, because that is what R does: it finishes
            // X_fixed %*% B and only then adds the Z product. Accumulating the p
            // terms straight into `column` would reassociate the sum.
            const double* b_gene = b.data + gene * p;
            for (std::int64_t i = 0; i < n; ++i) {
              double fixed = 0.0;
              for (std::int64_t k = 0; k < p; ++k) {
                fixed += x_fixed.data[i + k * n] * b_gene[k];
              }
              column[i] += fixed;
            }
          }
        }
      },
      interrupted);
  return status;
}

}  // namespace pace
