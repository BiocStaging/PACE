// gene_solve.hpp -- the penalised weighted least squares solve, one system per
// gene, for a chunk of genes.
//
// This is the numerical core of the IRLS iteration. For gene g it solves
//
//   [ X'WX   X'WZ ] [ beta ]   [ X'Wz ]
//   [ Z'WX   Z'WZ + Lambda ] [ u ] = [ Z'Wz ]
//
// with W = diag(w[, g]), Lambda = diag(lam[, g]) the 1/tau ridge, and returns
// the diagonal of the inverse of the whole system (the posterior variances the
// standard errors and the EM update read).
//
// Z is never formed. Each random-effect block assigns every cell to exactly one
// group, so Z'WZ is block diagonal within a block, and the only dense coupling
// is between blocks; the work is done on the per-group term matrices instead.
// With two blocks the solve is Schur-partitioned on that structure, which is
// what makes the per-gene cost O(K_g K_t^3 + q1 q2^2 + p q^2) rather than
// O((p + q)^3).
//
// Determinism: the group loop and the gene loop each write only their own
// output slots, so the result does not depend on the thread count.
#ifndef PACE_GENE_SOLVE_HPP
#define PACE_GENE_SOLVE_HPP

#include <cstdint>
#include <vector>

#include "core_types.hpp"

namespace pace {

// One random-effect block: its terms, which group each cell belongs to, and the
// cells of each group. `terms` is n * n_terms column-major; `cell_group` is
// 1-based with n entries; `group_cells` holds the 0-based cell indices of group
// g in [group_offsets[g], group_offsets[g + 1]).
struct SolveBlock {
  int col_offset = 0;
  int n_terms = 0;
  int n_groups = 0;
  Span<const double> terms;
  Span<const int> cell_group;
  Span<const int> group_offsets;
  Span<const int> group_cells;
};

// Solve one chunk of genes.
//
// `single_precision` runs the interior in float, as the R caller's
// `interior_precision = 1` did: the same expressions in single precision, about
// 1.5x faster on this hardware, with the final iteration run in double so the
// standard errors are full precision. A gene whose system is not positive
// definite leaves its column non-finite, which the caller detects and re-solves.
//
// Shapes: `x_fixed` is n * p column-major; `w` and `z` are n * n_genes;
// `lam_diag` is q_total * n_genes; `beta_out` is p * n_genes; `u_out` is
// q_total * n_genes; `ainv_diag_out` is (p + q_total) * n_genes.
Status solve_genes_chunk(Span<const double> x_fixed, std::int64_t n, std::int64_t p,
                         const std::vector<SolveBlock>& blocks, Span<const double> w,
                         Span<const double> z, Span<const double> lam_diag,
                         std::int64_t q_total, std::int64_t n_genes, bool single_precision,
                         int n_threads, const InterruptCheck& interrupted, Span<double> beta_out,
                         Span<double> u_out, Span<double> ainv_diag_out);

}  // namespace pace

#endif  // PACE_GENE_SOLVE_HPP
