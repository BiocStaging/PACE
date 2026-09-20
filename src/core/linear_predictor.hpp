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
// caller's buffer, in parallel.
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
//
// ---------------------------------------------------------------------------
// Two ways of reading the same design
// ---------------------------------------------------------------------------
//
// Z is never a general sparse matrix. It is a stack of blocks, and in a block
// every cell belongs to exactly one group, so
//
//   Z[i, col_offset + t * n_groups + cell_group[i]] = terms[i, t]
//
// and every other entry of the block's rows is zero. Walking Z as a CSC
// therefore streams the whole design ONCE PER GENE -- 4 + 8 bytes per non-zero,
// 338 MB a gene at 1.2M cells and 23 non-zeros a row -- and touches the gene's
// n-long eta column once per (term, group), which at that size is far past the
// last level of cache.
//
// Given the block structure the same sum can be taken cell-panel by cell-panel
// with the genes on the INSIDE. The design panel is then read once per gene
// CHUNK instead of once per gene, and the eta panel stays in cache while the
// terms are accumulated into it. That is what `blocks` selects below. The two
// paths are bit-identical, not merely close: a cell's contributions arrive in
// the same order either way (blocks in column order, terms ascending within a
// block), and the fixed part is still added only once the Z sum is complete.
#ifndef PACE_LINEAR_PREDICTOR_HPP
#define PACE_LINEAR_PREDICTOR_HPP

#include <cstdint>
#include <vector>

#include "core_types.hpp"
#include "count_stats.hpp"
#include "gene_solve.hpp"

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
// `blocks` is the same design seen as the stack of blocks that built it (the
// solve's own view of it). Pass it whenever the caller has it: the result is
// identical to the bit and the design is read once per chunk rather than once
// per gene. Pass nullptr to walk `z` instead, which is all a caller without the
// block structure can do. When `blocks` is given, `z` is used only for its
// shape, and the blocks must tile Z's columns in order -- block b starting at
// `col_offset` and holding n_terms * n_groups of them.
//
// `cell_panel` is how many cells one panel of the blocked path covers. Zero
// picks a panel from the design's shape, which is what every caller in the
// package does; a positive value overrides it, which is how the panel size is
// measured. It is ignored when `blocks` is nullptr.
//
// Threading: the CSC path splits the genes and the blocked path splits the
// cells, and in both every output element is written by exactly one worker, so
// the result does not depend on `n_threads`.
Status eta_block(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                 std::int64_t p, Span<const double> b, const CscView& z,
                 const std::vector<SolveBlock>* blocks, Span<const double> u,
                 Span<const int> genes, std::int64_t n, std::int64_t n_genes_total,
                 Span<double> eta, std::int64_t cell_panel, int n_threads,
                 const InterruptCheck& interrupted);

}  // namespace pace

#endif  // PACE_LINEAR_PREDICTOR_HPP
