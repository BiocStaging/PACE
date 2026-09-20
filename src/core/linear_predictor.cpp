// fp_no_contract.hpp must come first: it forbids FMA contraction for the rest of
// the translation unit, which is what keeps this bit-identical to the R it
// replaces on FMA targets.
#include "fp_no_contract.hpp"

#include "linear_predictor.hpp"

#include <algorithm>
#include <atomic>
#include <cstring>

#include "thread_pool.hpp"

namespace pace {

namespace {

// The fixed-effect part of one gene, added over the half-open cell range
// [first, last) once that range's Z sum is complete.
//
// `b_gene` points at the gene's p coefficients. For p > 1 each element's fixed
// part is summed over p on its own before it reaches the accumulator, because
// that is what R does: it finishes X_fixed %*% B and only then adds the Z
// product. Accumulating the p terms straight into `column` would reassociate
// the sum.
void add_fixed_part(const double* x1, bool x1_is_unit, const double* x_fixed, std::int64_t p,
                    const double* b_gene, std::int64_t n, std::int64_t first, std::int64_t last,
                    double* column) {
  if (p == 1) {
    const double coefficient = b_gene[0];
    if (x1_is_unit) {
      for (std::int64_t i = first; i < last; ++i) column[i] += coefficient;
    } else {
      for (std::int64_t i = first; i < last; ++i) column[i] += x1[i] * coefficient;
    }
    return;
  }
  for (std::int64_t i = first; i < last; ++i) {
    double fixed = 0.0;
    for (std::int64_t k = 0; k < p; ++k) {
      fixed += x_fixed[i + k * n] * b_gene[k];
    }
    column[i] += fixed;
  }
}

// Cells per panel of the blocked path, when the caller has not chosen one.
//
// While a panel is being filled a worker touches four genes' eta panels, the
// panel's slice of every term column of every block, and the panel's slice of
// the fixed design. The panel wants to be small enough for that set to sit in
// the first two levels of cache and large enough that the loops are long.
//
// AD-HOC: measured over 256 to 8192 at three shapes (breast cancer, melanoma,
// and a 1.2M cell x 5,001 gene panel) and at gene chunks of 16 and 128. The
// curve is flat between 512 and 1024 everywhere, and above 1536 the wide-chunk
// cases fall away, so this sits at the top of the flat part. A multiple of 64
// doubles puts a panel boundary on a cache-line boundary whenever the buffer
// itself is aligned, which keeps two workers off the same line.
const std::int64_t kCellPanel = 1024;

// eta by walking Z as a compressed-column matrix: the only thing a caller
// without the block structure can do. Genes are split across workers.
Status eta_block_csc(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                     std::int64_t p, Span<const double> b, const CscView& z,
                     Span<const double> u, Span<const int> genes, std::int64_t n,
                     Span<double> eta, int n_threads, const InterruptCheck& interrupted) {
  const std::int64_t n_chunk = genes.size;
  const std::int64_t q = z.n_cols;
  return parallel_for(
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
          add_fixed_part(x1.data, x1_is_unit, x_fixed.data, p, b.data + gene * p, n, 0, n, column);
        }
      },
      interrupted);
}

// eta from the block structure: a panel of cells at a time, genes on the inside.
//
// Cell i of block b sits in group cell_group[i] - 1, so its contribution from
// column (t, g) of that block is terms[i, t] * u[col_offset + t * n_groups + g]
// and it is the ONLY column of the block that reaches it. Running the blocks in
// order and the terms ascending within a block therefore hands each cell its
// contributions in exactly the order the CSC walk does, because that is the
// order the columns were laid out in (random_design.cpp, and the cbind in
// build_random_design_multi).
//
// The CSC walk skips a column whose coefficient is exactly zero; this does not,
// and adds 0 * terms[i, t] instead. For finite terms that is bit-identical:
// x + 0 is x for every x the accumulator can hold here, and the accumulator
// starts at +0 so it is never -0. A non-finite entry in the design would break
// that, and would already have made the gene's solve non-finite.
//
// Cells are split across workers, so every eta element is still written by
// exactly one worker.
Status eta_block_blocked(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                         std::int64_t p, Span<const double> b,
                         const std::vector<SolveBlock>& blocks, std::int64_t q,
                         Span<const double> u, Span<const int> genes, std::int64_t n,
                         Span<double> eta, std::int64_t cell_panel, int n_threads,
                         const InterruptCheck& interrupted) {
  const std::int64_t n_chunk = genes.size;
  const std::int64_t panel = cell_panel > 0 ? cell_panel : kCellPanel;
  const std::int64_t n_panels = (n + panel - 1) / panel;
  std::atomic<bool> group_out_of_range{false};

  const Status status = parallel_for(
      n_panels, n_threads, 1,
      [&](std::int64_t first_panel, std::int64_t last_panel) {
        for (std::int64_t panel_index = first_panel; panel_index < last_panel; ++panel_index) {
          const std::int64_t lo = panel_index * panel;
          const std::int64_t hi = std::min(n, lo + panel);

          // Every group index this panel will use, checked once rather than
          // once per gene. A bad one would index outside the gene's u column.
          for (const SolveBlock& block : blocks) {
            for (std::int64_t i = lo; i < hi; ++i) {
              const int group = block.cell_group[i] - 1;
              if (group < 0 || group >= block.n_groups) {
                group_out_of_range.store(true);
                return;
              }
            }
          }

          // ---- the Z part, in Z's own column order -----------------------
          const std::size_t panel_bytes = static_cast<std::size_t>(hi - lo) * sizeof(double);
          for (std::int64_t j = 0; j < n_chunk; ++j) {
            std::memset(eta.data + j * n + lo, 0, panel_bytes);
          }
          for (const SolveBlock& block : blocks) {
            const int* cell_group = block.cell_group.data;
            for (std::int64_t t = 0; t < block.n_terms; ++t) {
              const double* term = block.terms.data + t * n;
              const std::int64_t column_of_term = block.col_offset + t * block.n_groups;
              // Four genes at a time. term[i] and cell_group[i] are then loaded
              // once for the four rather than once each, which is what keeps
              // this loop issuing multiply-adds instead of waiting on loads.
              // Genes do not interact, so this is the same arithmetic in the
              // same order as one gene at a time.
              std::int64_t j = 0;
              for (; j + 4 <= n_chunk; j += 4) {
                double* column_0 = eta.data + j * n;
                double* column_1 = column_0 + n;
                double* column_2 = column_0 + 2 * n;
                double* column_3 = column_0 + 3 * n;
                const double* coefficient_0 = u.data + genes[j] * q + column_of_term;
                const double* coefficient_1 = u.data + genes[j + 1] * q + column_of_term;
                const double* coefficient_2 = u.data + genes[j + 2] * q + column_of_term;
                const double* coefficient_3 = u.data + genes[j + 3] * q + column_of_term;
                for (std::int64_t i = lo; i < hi; ++i) {
                  const double value = term[i];
                  const int group = cell_group[i] - 1;
                  column_0[i] += value * coefficient_0[group];
                  column_1[i] += value * coefficient_1[group];
                  column_2[i] += value * coefficient_2[group];
                  column_3[i] += value * coefficient_3[group];
                }
              }
              for (; j < n_chunk; ++j) {
                double* column = eta.data + j * n;
                const double* coefficient = u.data + genes[j] * q + column_of_term;
                for (std::int64_t i = lo; i < hi; ++i) {
                  column[i] += term[i] * coefficient[cell_group[i] - 1];
                }
              }
            }
          }

          // ---- the fixed part, added once the Z sum is complete -----------
          for (std::int64_t j = 0; j < n_chunk; ++j) {
            add_fixed_part(x1.data, x1_is_unit, x_fixed.data, p, b.data + genes[j] * p, n, lo, hi,
                           eta.data + j * n);
          }
        }
      },
      interrupted);
  if (!status.is_ok()) return status;
  if (group_out_of_range.load()) {
    return Status::failure(StatusCode::invalid_argument,
                           "eta_block: a cell's group is outside its block.");
  }
  return Status::success();
}

}  // namespace

Status eta_block(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                 std::int64_t p, Span<const double> b, const CscView& z,
                 const std::vector<SolveBlock>* blocks, Span<const double> u,
                 Span<const int> genes, std::int64_t n, std::int64_t n_genes_total,
                 Span<double> eta, std::int64_t cell_panel, int n_threads,
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

  if (blocks == nullptr) {
    return eta_block_csc(x1, x1_is_unit, x_fixed, p, b, z, u, genes, n, eta, n_threads,
                         interrupted);
  }

  // The blocks must be the columns of THIS Z, in this order: the fast path
  // reads the design through them and never looks at Z again.
  std::int64_t offset = 0;
  for (const SolveBlock& block : *blocks) {
    if (block.n_terms < 0 || block.n_groups < 1)
      return Status::failure(StatusCode::invalid_argument,
                             "eta_block: a block has no terms or no groups.");
    if (block.col_offset != offset)
      return Status::failure(StatusCode::invalid_argument,
                             "eta_block: the blocks do not tile Z's columns in order.");
    if (block.terms.size != n * block.n_terms)
      return Status::failure(StatusCode::invalid_argument,
                             "eta_block: a block's terms are not n * n_terms.");
    if (block.cell_group.size != n)
      return Status::failure(StatusCode::invalid_argument,
                             "eta_block: a block's cell_group does not have one entry per cell.");
    offset += static_cast<std::int64_t>(block.n_terms) * block.n_groups;
  }
  if (offset != q)
    return Status::failure(StatusCode::invalid_argument,
                           "eta_block: the blocks do not cover Z's columns.");

  return eta_block_blocked(x1, x1_is_unit, x_fixed, p, b, *blocks, q, u, genes, n, eta, cell_panel,
                           n_threads, interrupted);
}

}  // namespace pace
