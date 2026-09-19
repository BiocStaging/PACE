#ifndef PACE_CORE_RANDOM_DESIGN_HPP
#define PACE_CORE_RANDOM_DESIGN_HPP

#include <cstdint>

#include "core_types.hpp"

namespace pace {

// The sparse random-effect design Z of one block, built straight into
// compressed-column form.
//
// A block crosses K_terms model-matrix columns with K_groups groups, so its
// column (t, g) holds, for every cell i in group g, the value X_terms[i, t]:
//
//   Z[i, t * K_groups + g] = X_terms[i, t]   when cell_group[i] == g
//
// Every cell belongs to exactly one group, so each cell contributes exactly
// K_terms entries and the block has n * K_terms of them. Explicit zeros are
// KEPT, because Matrix::sparseMatrix() kept them and the solver's column
// pointers are what they are because of it.
//
// R built this from triplets: three vectors of n * K_terms elements each --
// at 1.2M cells and twenty terms roughly 400 MB of allocation -- handed to
// Matrix::sparseMatrix(), which then SORTED them into column order. None of
// that is necessary. Walking groups in order and cells within a group in
// ascending order emits the rows already sorted, so this writes the final CSC
// arrays once, with no triplet buffers and no sort.
//
// `cell_group` holds a 0-based group per cell. `column_pointer` must have
// K_terms * K_groups + 1 entries, and `row_index` and `values` must each have
// n * K_terms. `cells_by_group` and `group_start` are the cell indices grouped
// by group, ascending within a group, which the caller also wants.
Status random_design_block(Span<const double> x_terms, Span<const int> cell_group,
                           std::int64_t n, std::int64_t n_terms, std::int64_t n_groups,
                           Span<int> column_pointer, Span<int> row_index, Span<double> values,
                           Span<int> cells_by_group, Span<int> group_start);

}  // namespace pace

#endif  // PACE_CORE_RANDOM_DESIGN_HPP
