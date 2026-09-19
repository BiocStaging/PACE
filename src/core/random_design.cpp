#include "random_design.hpp"

#include <algorithm>
#include <vector>

namespace pace {

Status random_design_block(Span<const double> x_terms, Span<const int> cell_group,
                           std::int64_t n, std::int64_t n_terms, std::int64_t n_groups,
                           Span<int> column_pointer, Span<int> row_index, Span<double> values,
                           Span<int> cells_by_group, Span<int> group_start) {
  const std::int64_t n_columns = n_terms * n_groups;
  const std::int64_t nnz = n * n_terms;
  if (x_terms.size != n * n_terms) {
    return Status::failure(StatusCode::invalid_argument, "x_terms must be n * n_terms");
  }
  if (cell_group.size != n) {
    return Status::failure(StatusCode::invalid_argument, "cell_group must have one entry per cell");
  }
  if (column_pointer.size != n_columns + 1 || row_index.size != nnz || values.size != nnz) {
    return Status::failure(StatusCode::invalid_argument, "the output arrays have the wrong size");
  }
  if (cells_by_group.size != n || group_start.size != n_groups + 1) {
    return Status::failure(StatusCode::invalid_argument, "the grouping arrays have the wrong size");
  }

  // Count each group, then turn the counts into starts. A group with no cells
  // is allowed: it contributes an empty column, which is what split() gave the
  // R code as integer(0), and the solver relies on the column still existing.
  for (std::int64_t g = 0; g <= n_groups; ++g) group_start[g] = 0;
  for (std::int64_t i = 0; i < n; ++i) {
    const int g = cell_group[i];
    if (g < 0 || g >= n_groups) {
      return Status::failure(StatusCode::invalid_argument, "a cell's group is out of range");
    }
    group_start[g + 1] += 1;
  }
  for (std::int64_t g = 0; g < n_groups; ++g) group_start[g + 1] += group_start[g];

  // Fill the group membership in ONE pass over the cells in ascending order,
  // which is what leaves each group's cells sorted and lets the columns below
  // be written without a sort.
  {
    std::vector<int> fill_position(static_cast<std::size_t>(n_groups));
    for (std::int64_t g = 0; g < n_groups; ++g) {
      fill_position[static_cast<std::size_t>(g)] = group_start[g];
    }
    for (std::int64_t i = 0; i < n; ++i) {
      const int g = cell_group[i];
      cells_by_group[fill_position[static_cast<std::size_t>(g)]++] = static_cast<int>(i);
    }
  }

  // Column (t, g) is at t * n_groups + g and holds that group's cells, so its
  // length is the group's size whatever the term.
  column_pointer[0] = 0;
  for (std::int64_t t = 0; t < n_terms; ++t) {
    for (std::int64_t g = 0; g < n_groups; ++g) {
      const std::int64_t column = t * n_groups + g;
      const std::int64_t group_size = group_start[g + 1] - group_start[g];
      column_pointer[column + 1] = column_pointer[column] + static_cast<int>(group_size);
    }
  }

  for (std::int64_t t = 0; t < n_terms; ++t) {
    const double* term_column = x_terms.data + t * n;
    for (std::int64_t g = 0; g < n_groups; ++g) {
      const std::int64_t column = t * n_groups + g;
      std::int64_t position = column_pointer[column];
      for (std::int64_t k = group_start[g]; k < group_start[g + 1]; ++k) {
        const int cell = cells_by_group[k];
        row_index[position] = cell;
        values[position] = term_column[cell];
        ++position;
      }
    }
  }
  return Status::success();
}

}  // namespace pace
