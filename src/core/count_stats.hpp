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
// belong to no group.
//
// Shapes: group has n_rows entries with codes in [-1, n_groups). `means` is
// caller-owned with n_groups * n_cols entries, column-major (groups x genes).
Status group_column_means(const CscView& counts, Span<const int> group, int n_groups,
                          bool detection, Span<double> means, int n_threads,
                          const InterruptCheck& interrupted);

}  // namespace pace

#endif  // PACE_COUNT_STATS_HPP
