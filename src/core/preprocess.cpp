// preprocess.cpp -- implementation of preprocess.hpp.
#include "fp_no_contract.hpp"  // must precede the arithmetic below

#include "preprocess.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include "r_summaries.hpp"

namespace pace {

Status drop_sparse_kernel(Span<double> kernel, std::int64_t n, std::int64_t n_types,
                          Span<const int> celltype, double min_effective,
                          std::int64_t* n_dropped) {
  if (kernel.size != n * n_types || celltype.size != n) {
    return Status::failure(StatusCode::invalid_argument, "kernel must be n x n_types");
  }
  std::vector<std::vector<std::int64_t>> cells_of_type(static_cast<std::size_t>(n_types));
  for (std::int64_t i = 0; i < n; ++i) {
    const int type = celltype[i];
    if (type >= 0 && type < n_types) cells_of_type[static_cast<std::size_t>(type)].push_back(i);
  }
  std::int64_t dropped = 0;
  for (std::int64_t focal = 0; focal < n_types; ++focal) {
    const std::vector<std::int64_t>& cells = cells_of_type[static_cast<std::size_t>(focal)];
    if (cells.empty()) continue;
    for (std::int64_t neighbour = 0; neighbour < n_types; ++neighbour) {
      double* column = kernel.data + neighbour * n;
      const double mean = r_mean(column, cells, true, nullptr);
      long double total = 0.0L;
      double largest = -std::numeric_limits<double>::infinity();
      for (std::int64_t cell : cells) {
        const double centred = column[cell] - mean;
        const double square = centred * centred;
        if (std::isnan(square)) continue;
        total += square;
        if (square > largest) largest = square;
      }
      const double effective = largest > 0 ? static_cast<double>(total) / largest : 0.0;
      if (effective < min_effective) {
        for (std::int64_t cell : cells) column[cell] = 0.0;
        dropped += 1;
      }
    }
  }
  *n_dropped = dropped;
  return Status::success();
}

Status centre_within_groups(Span<double> kernel, std::int64_t n, std::int64_t n_types,
                            Span<const int> image, Span<const int> celltype, int n_images,
                            int n_celltypes) {
  if (kernel.size != n * n_types || image.size != n || celltype.size != n) {
    return Status::failure(StatusCode::invalid_argument, "kernel must be n x n_types");
  }
  // One group per (image, cell type) pair, as interaction() formed them.
  const std::int64_t n_groups = static_cast<std::int64_t>(n_images) * n_celltypes;
  std::vector<std::vector<std::int64_t>> cells_of_group(static_cast<std::size_t>(n_groups) + 1);
  for (std::int64_t i = 0; i < n; ++i) {
    const int im = image[i];
    const int type = celltype[i];
    const std::int64_t group =
        (im >= 0 && im < n_images && type >= 0 && type < n_celltypes)
            ? static_cast<std::int64_t>(im) + static_cast<std::int64_t>(type) * n_images
            : n_groups;
    cells_of_group[static_cast<std::size_t>(group)].push_back(i);
  }
  for (std::int64_t column_index = 0; column_index < n_types; ++column_index) {
    double* column = kernel.data + column_index * n;
    for (const std::vector<std::int64_t>& cells : cells_of_group) {
      if (cells.empty()) continue;
      const double mean = r_mean(column, cells, false, nullptr);
      for (std::int64_t cell : cells) column[cell] = column[cell] - mean;
    }
  }
  return Status::success();
}

Status standardise_column(Span<double> values, std::int64_t n, double* standard_deviation) {
  if (values.size != n) {
    return Status::failure(StatusCode::invalid_argument, "values must have n entries");
  }
  if (n == 0) return Status::success();
  // scale() centres with colMeans(), which is a long double sum over n.
  long double sum = 0.0L;
  for (std::int64_t i = 0; i < n; ++i) sum += values[i];
  const double centre = static_cast<double>(sum / static_cast<long double>(n));
  long double squares = 0.0L;
  for (std::int64_t i = 0; i < n; ++i) {
    values[i] = values[i] - centre;
    const double square = values[i] * values[i];
    squares += square;
  }
  const std::int64_t denominator = n - 1 > 1 ? n - 1 : 1;
  const double spread =
      std::sqrt(static_cast<double>(squares / static_cast<long double>(denominator)));
  for (std::int64_t i = 0; i < n; ++i) values[i] = values[i] / spread;
  if (standard_deviation != nullptr) *standard_deviation = spread;
  return Status::success();
}

Status anchor_mask(Span<const double> core_means, std::int64_t n_types, std::int64_t n_genes,
                   double owner_threshold, double core_threshold, Span<double> mask,
                   Span<int> n_anchor) {
  if (core_means.size != n_types * n_genes || mask.size != n_types * n_genes ||
      n_anchor.size != n_types) {
    return Status::failure(StatusCode::invalid_argument, "core_means must be n_types x n_genes");
  }
  for (std::int64_t t = 0; t < n_types; ++t) n_anchor[t] = 0;
  for (std::int64_t gene = 0; gene < n_genes; ++gene) {
    const double* column = core_means.data + gene * n_types;
    double owner_mean = column[0];
    std::int64_t owner = 0;
    for (std::int64_t t = 1; t < n_types; ++t) {
      if (column[t] > owner_mean) {  // the first maximum wins, as which.max does
        owner_mean = column[t];
        owner = t;
      }
    }
    for (std::int64_t t = 0; t < n_types; ++t) {
      const bool is_anchor = owner != t && owner_mean > owner_threshold &&
                             column[t] / r_pmax(owner_mean, 1e-9) < core_threshold;
      mask[t + gene * n_types] = is_anchor ? 1.0 : 0.0;
      if (is_anchor) n_anchor[t] += 1;
    }
  }
  return Status::success();
}

Status normalise_rows_to_max(Span<double> values, std::int64_t n_rows, std::int64_t n_cols) {
  if (values.size != n_rows * n_cols) {
    return Status::failure(StatusCode::invalid_argument, "values must be n_rows x n_cols");
  }
  for (std::int64_t row = 0; row < n_rows; ++row) {
    double largest = -std::numeric_limits<double>::infinity();
    for (std::int64_t col = 0; col < n_cols; ++col) {
      const double value = values[row + col * n_rows];
      if (std::isnan(value)) continue;
      if (value > largest) largest = value;
    }
    const bool usable = std::isfinite(largest) && largest > 0;
    for (std::int64_t col = 0; col < n_cols; ++col) {
      double value = values[row + col * n_rows];
      value = usable ? value / largest : 0.0;
      if (std::isnan(value)) value = 0.0;
      values[row + col * n_rows] = value > 1.0 ? 1.0 : value;
    }
  }
  return Status::success();
}

}  // namespace pace
