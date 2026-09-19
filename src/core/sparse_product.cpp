#include "sparse_product.hpp"

#include <algorithm>
#include <cstring>

#include "thread_pool.hpp"

namespace pace {
namespace {

// One column of the product, accumulated densely and read back through the list
// of rows that were actually touched. `marker` records which rows are live for
// this column without clearing the whole scratch between columns.
struct SparseAccumulator {
  std::vector<double> scratch;
  std::vector<int> touched;
  std::vector<char> marker;

  void resize(std::int64_t n_rows) {
    if (static_cast<std::int64_t>(scratch.size()) == n_rows) return;
    scratch.assign(static_cast<std::size_t>(n_rows), 0.0);
    marker.assign(static_cast<std::size_t>(n_rows), 0);
    touched.reserve(static_cast<std::size_t>(n_rows));
  }

  // W times one column of the right operand, left in `scratch` with the live
  // rows listed in `touched`, sorted so the output column is in row order.
  void accumulate(const CscView& left, const CscView& right, std::int64_t column) {
    touched.clear();
    for (int p = right.column_pointer[column]; p < right.column_pointer[column + 1]; ++p) {
      const int right_row = right.row_index[p];
      const double right_value = right.values[p];
      if (right_value == 0.0) continue;
      for (int q = left.column_pointer[right_row]; q < left.column_pointer[right_row + 1]; ++q) {
        const int row = left.row_index[q];
        if (!marker[static_cast<std::size_t>(row)]) {
          marker[static_cast<std::size_t>(row)] = 1;
          touched.push_back(row);
        }
        scratch[static_cast<std::size_t>(row)] += left.values[q] * right_value;
      }
    }
    std::sort(touched.begin(), touched.end());
  }

  void clear() {
    for (std::size_t k = 0; k < touched.size(); ++k) {
      scratch[static_cast<std::size_t>(touched[k])] = 0.0;
      marker[static_cast<std::size_t>(touched[k])] = 0;
    }
    touched.clear();
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
      accumulator.clear();
    }
  };
  Status status = parallel_for(n_columns, n_threads, 1, count_body, interrupted);
  if (!status.is_ok()) return status;

  column_pointer.assign(static_cast<std::size_t>(n_columns) + 1, 0);
  for (std::int64_t j = 0; j < n_columns; ++j) {
    column_pointer[static_cast<std::size_t>(j) + 1] =
        column_pointer[static_cast<std::size_t>(j)] + column_count[static_cast<std::size_t>(j)];
  }
  const std::int64_t total = column_pointer[static_cast<std::size_t>(n_columns)];
  row_index.assign(static_cast<std::size_t>(total), 0);
  values.assign(static_cast<std::size_t>(total), 0.0);

  // Pass two: fill. Column j writes only into its own slice, so this is the
  // same partition as the count pass and the result does not depend on it.
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
      accumulator.clear();
    }
  };
  return parallel_for(n_columns, n_threads, 1, fill_body, interrupted);
}

}  // namespace pace
