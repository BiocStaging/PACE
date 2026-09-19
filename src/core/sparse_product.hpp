#ifndef PACE_CORE_SPARSE_PRODUCT_HPP
#define PACE_CORE_SPARSE_PRODUCT_HPP

#include <cstdint>
#include <vector>

#include "core_types.hpp"
#include "count_stats.hpp"  // CscView

namespace pace {

// The ambient cache A = W Y, where W is the n x n cross-cell-type weight matrix
// and Y the n x G counts, both compressed-column.
//
// Column j of A is one sparse mat-vec, W times column j of Y, and the columns
// are independent, so this parallelises over genes and each worker owns its own
// accumulator. The accumulator is Gustavson's: a dense scratch column plus a
// list of the rows touched, which makes each column linear in its own non-zeros
// instead of in n.
//
// Structural zeros are DROPPED, matching what Matrix's %*% hands back -- an
// entry that cancels exactly is not stored, and the solver reads this one gene
// at a time, so a stored zero would only cost it a wasted iteration.
//
// Two passes: the first counts each column's non-zeros so the column pointers
// can be laid out before anything is written, the second fills. That keeps the
// output exactly the size it needs, which matters because R's version
// materialised the whole product and then converted it -- two full n x G
// allocations for a matrix the core only ever reads column by column.
// `first_column` and `n_columns` select a CONTIGUOUS RANGE of the right
// operand's columns, which is how the solver consumes it: one gene chunk at a
// time. The output has exactly `n_columns` columns, numbered from zero, so it
// drops straight into a GeneBlock with first_gene = 0.
//
// The range matters for a reason beyond convenience. A sparse product is
// COLUMN-INDEPENDENT -- output column j depends only on input column j -- so
// computing a range gives bit-for-bit what computing everything and slicing
// would. Nothing about the answer depends on the chunking.
Status sparse_product_csc(const CscView& left, const CscView& right,
                          std::int64_t first_column, std::int64_t n_columns,
                          std::vector<int>& column_pointer, std::vector<int>& row_index,
                          std::vector<double>& values, int n_threads,
                          const InterruptCheck& interrupted);

}  // namespace pace

#endif  // PACE_CORE_SPARSE_PRODUCT_HPP
