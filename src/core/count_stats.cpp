// count_stats.cpp -- implementation of count_stats.hpp.
#include "count_stats.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "thread_pool.hpp"

namespace pace {

Status group_column_means(const CscView& counts, Span<const int> group, int n_groups,
                          bool detection, Span<double> means, int n_threads,
                          const InterruptCheck& interrupted) {
  if (group.size != counts.n_rows) {
    return Status::failure(StatusCode::invalid_argument, "group length differs from the number of cells");
  }
  if (counts.column_pointer.size != counts.n_cols + 1) {
    return Status::failure(StatusCode::invalid_argument, "column pointer length differs from n_cols + 1");
  }
  if (means.size != static_cast<std::int64_t>(n_groups) * counts.n_cols) {
    return Status::failure(StatusCode::invalid_argument, "output size differs from n_groups * n_cols");
  }

  for (std::int64_t k = 0; k < counts.values.size; ++k) {
    if (!std::isfinite(counts.values[k])) {
      return Status::failure(StatusCode::invalid_argument, "counts must be finite (no NA, NaN or Inf)");
    }
  }

  std::vector<std::int64_t> group_size(n_groups, 0);
  for (std::int64_t i = 0; i < group.size; ++i) {
    const int g = group[i];
    if (g >= 0 && g < n_groups) group_size[g] += 1;
  }

  auto body = [&](std::int64_t begin, std::int64_t end) {
    std::vector<long double> column_sum(n_groups);
    for (std::int64_t gene = begin; gene < end; ++gene) {
      std::fill(column_sum.begin(), column_sum.end(), 0.0L);
      for (int k = counts.column_pointer[gene]; k < counts.column_pointer[gene + 1]; ++k) {
        const int g = group[counts.row_index[k]];
        if (g < 0 || g >= n_groups) continue;
        const double value = counts.values[k];
        if (detection) {
          column_sum[g] += value > 0.0 ? 1.0L : 0.0L;
        } else {
          column_sum[g] += value;
        }
      }
      for (int g = 0; g < n_groups; ++g) {
        const long double mean = column_sum[g] / static_cast<long double>(group_size[g]);
        means[g + gene * static_cast<std::int64_t>(n_groups)] = static_cast<double>(mean);
      }
    }
  };
  return parallel_for(counts.n_cols, n_threads, 16, body, interrupted);
}

}  // namespace pace
