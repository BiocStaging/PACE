#include "fp_no_contract.hpp"  // must precede the arithmetic below

#include "sparse_product.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <limits>
#include <vector>

#include "thread_pool.hpp"

namespace pace {
namespace {

// One column of the product, accumulated densely and read back through the list
// of rows that were actually touched. `marker` records which rows are live for
// this column without clearing the whole scratch between columns.
struct SparseAccumulator {
  std::vector<double> scratch;
  std::vector<int> touched;
  // A GENERATION per row rather than a flag. A row belongs to the column being
  // built when its generation matches the current one, so a new column costs an
  // increment instead of clearing n entries, and `scratch` never has to be
  // zeroed at all: the first touch in a column assigns rather than accumulates.
  // At 1.2M cells that is 9.6 MB of memset saved per column.
  std::vector<std::uint32_t> generation;
  std::uint32_t current_generation = 0;

  void resize(std::int64_t n_rows) {
    if (static_cast<std::int64_t>(scratch.size()) == n_rows) return;
    scratch.resize(static_cast<std::size_t>(n_rows));
    generation.assign(static_cast<std::size_t>(n_rows), 0);
    current_generation = 0;
    // No reserve on `touched`: a chunk's live set is a small fraction of n, and
    // reserving n would allocate 4.8 MB per worker per call for nothing.
  }

  // W times one column of the right operand, left in `scratch` with the live
  // rows listed in `touched`, sorted so the output column is in row order.
  void accumulate(const CscView& left, const CscView& right, std::int64_t column) {
    touched.clear();
    ++current_generation;
    if (current_generation == 0) {          // wrapped: retire every stale mark
      std::fill(generation.begin(), generation.end(), 0);
      current_generation = 1;
    }
    for (int p = right.column_pointer[column]; p < right.column_pointer[column + 1]; ++p) {
      const int right_row = right.row_index[p];
      const double right_value = right.values[p];
      if (right_value == 0.0) continue;
      for (int q = left.column_pointer[right_row]; q < left.column_pointer[right_row + 1]; ++q) {
        const int row = left.row_index[q];
        const std::size_t slot = static_cast<std::size_t>(row);
        const double contribution = left.values[q] * right_value;
        if (generation[slot] != current_generation) {
          generation[slot] = current_generation;
          scratch[slot] = contribution;      // first touch: assign, do not add
          touched.push_back(row);
        } else {
          scratch[slot] += contribution;
        }
      }
    }
    std::sort(touched.begin(), touched.end());
  }
};

}  // namespace

Status sparse_product_csc(const CscView& left, const CscView& right,
                          std::int64_t first_column, std::int64_t n_columns,
                          std::vector<int>& column_pointer, std::vector<int>& row_index,
                          std::vector<double>& values, int n_threads,
                          const InterruptCheck& interrupted) {
  if (left.n_cols != right.n_rows) {
    return Status::failure(StatusCode::invalid_argument,
                           "the ambient field and the counts do not conform");
  }
  if (first_column < 0 || n_columns < 0 || first_column + n_columns > right.n_cols) {
    return Status::failure(StatusCode::invalid_argument, "the column range is out of bounds");
  }
  const std::int64_t n_rows = left.n_rows;

  // Pass one: how many non-zeros each column of the product holds. Counting
  // first means the column pointers are known before a single value is written,
  // so the output is allocated once at exactly its final size.
  std::vector<int> column_count(static_cast<std::size_t>(n_columns), 0);
  auto count_body = [&](std::int64_t begin, std::int64_t end) {
    static thread_local SparseAccumulator accumulator;
    accumulator.resize(n_rows);
    for (std::int64_t j = begin; j < end; ++j) {
      accumulator.accumulate(left, right, first_column + j);
      int live = 0;
      for (std::size_t k = 0; k < accumulator.touched.size(); ++k) {
        if (accumulator.scratch[static_cast<std::size_t>(accumulator.touched[k])] != 0.0) ++live;
      }
      column_count[static_cast<std::size_t>(j)] = live;
    }
  };
  Status status = parallel_for(n_columns, n_threads, 1, count_body, interrupted);
  if (!status.is_ok()) return status;

  // The prefix sum accumulates in int64 and is checked before it is written
  // back into the int column pointers. A compressed-column matrix cannot hold
  // more than INT_MAX non-zeros, and overflowing silently would give a negative
  // total, a wrapped allocation and a heap overrun rather than an error. At
  // 1.2M cells this needs a chunk of roughly 1,800 fully dense columns, so it
  // is far from the default 128 -- but nothing else checks it.
  column_pointer.assign(static_cast<std::size_t>(n_columns) + 1, 0);
  std::int64_t running_total = 0;
  for (std::int64_t j = 0; j < n_columns; ++j) {
    running_total += column_count[static_cast<std::size_t>(j)];
    if (running_total > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
      return Status::failure(StatusCode::invalid_argument,
                             "the chunk's ambient product exceeds what a sparse matrix can hold; "
                             "reduce chunk_size");
    }
    column_pointer[static_cast<std::size_t>(j) + 1] = static_cast<int>(running_total);
  }
  const std::int64_t total = running_total;
  // resize, not assign: pass two writes exactly `total` entries, so
  // value-initialising them first is pure waste -- at 1.2M cells and a 30%
  // dense chunk that is half a gigabyte of memset per product, and there are
  // hundreds of products in a fit.
  row_index.resize(static_cast<std::size_t>(total));
  values.resize(static_cast<std::size_t>(total));

  // Pass two: fill. Column j writes only into its own slice, so this is the
  // same partition as the count pass and the result does not depend on it.
  std::atomic<bool> failed(false);
  auto fill_body = [&](std::int64_t begin, std::int64_t end) {
    static thread_local SparseAccumulator accumulator;
    accumulator.resize(n_rows);
    for (std::int64_t j = begin; j < end; ++j) {
      accumulator.accumulate(left, right, first_column + j);
      std::int64_t position = column_pointer[static_cast<std::size_t>(j)];
      for (std::size_t k = 0; k < accumulator.touched.size(); ++k) {
        const int row = accumulator.touched[k];
        const double value = accumulator.scratch[static_cast<std::size_t>(row)];
        if (value == 0.0) continue;
        row_index[static_cast<std::size_t>(position)] = row;
        values[static_cast<std::size_t>(position)] = value;
        ++position;
      }
      // The two passes run the same accumulation, so they cannot disagree for a
      // non-negative W and Y -- no cancellation is reachable. This is written as
      // a general utility though, and a caller with signed weights would get a
      // silent write past the column's end. Cheap to turn that into an error.
      if (position != column_pointer[static_cast<std::size_t>(j) + 1]) {
        failed.store(true);
      }
    }
  };
  const Status fill_status = parallel_for(n_columns, n_threads, 1, fill_body, interrupted);
  if (!fill_status.is_ok()) return fill_status;
  if (failed.load()) {
    return Status::failure(StatusCode::internal_error,
                           "the counting and filling passes of the sparse product disagree");
  }
  return Status::success();
}

}  // namespace pace
