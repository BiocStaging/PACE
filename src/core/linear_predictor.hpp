// linear_predictor.hpp -- eta for one gene chunk.
//
// This replaces the R closure `.eta_block()`:
//
//   .xb_chunk(B, genes) + as.matrix(Z %*% U[, genes, drop = FALSE])
//
// which was the last arithmetic left in the streaming solver's hot loop. In R
// it allocated three dense n x chunk matrices per call (the sparse product, the
// as.matrix() copy, and the sum) -- 626 MB each at 1.2M cells and a chunk of
// 64 -- and ran single-threaded. Here the result is written once, into the
// caller's buffer, in parallel over genes.
//
// Summation order follows the R it replaces: the Z contributions are
// accumulated first, in column order of Z exactly as a CSC-times-dense product
// does, and the fixed-effect part is added afterwards. (Elementwise addition of
// the two finished parts is commutative, so which one R adds first does not
// matter; the order *within* each part does.)
//
// Agreement, measured over 432 shapes against the R expression:
//   p == 1  bit-identical -- the fixed part is one broadcast multiply
//   p  > 1  ~2e-15        -- R sends X_fixed %*% B through BLAS, whose
//                            accumulation and vectorisation are not reproducible
//                            portably; this sums each element over p in order.
// Both are far inside the 1e-10 the fixtures are gated at.
#ifndef PACE_LINEAR_PREDICTOR_HPP
#define PACE_LINEAR_PREDICTOR_HPP

#include <cstdint>

#include "core_types.hpp"
#include "count_stats.hpp"

namespace pace {

// eta[, j] = X_fixed %*% b[, genes[j]] + Z %*% u[, genes[j]]
//
// Fixed part: when `p == 1` the single column is `x1` (length n) and the
// contribution is x1[i] * b[genes[j]], with `x1_is_unit` skipping the multiply;
// `x_fixed` is then unused and may be empty. When `p > 1`, `x_fixed` is the
// dense n * p column-major design and `x1` is unused.
//
// Shapes: `b` is p * n_genes_total column-major, `u` is q * n_genes_total
// column-major, both indexed by the ORIGINAL gene number in `genes` (0-based).
// `z` is the n x q random-effect design. `eta` is n * n_chunk column-major and
// is fully overwritten.
//
// Threading: gene j is written only by the worker that owns it, so the result
// does not depend on `n_threads`.
Status eta_block(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                 std::int64_t p, Span<const double> b, const CscView& z,
                 Span<const double> u, Span<const int> genes, std::int64_t n,
                 std::int64_t n_genes_total, Span<double> eta, int n_threads,
                 const InterruptCheck& interrupted);

}  // namespace pace

#endif  // PACE_LINEAR_PREDICTOR_HPP
