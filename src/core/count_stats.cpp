// count_stats.cpp -- implementation of count_stats.hpp.
#include "count_stats.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
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

Status group_covariances(Span<const double> values, std::int64_t n, std::int64_t p,
                         Span<const int> group, int n_groups, Span<double> covariance,
                         int n_threads, const InterruptCheck& interrupted) {
  if (values.size != n * p || group.size != n) {
    return Status::failure(StatusCode::invalid_argument, "values must be n * p and group must have n entries");
  }
  if (covariance.size != static_cast<std::int64_t>(n_groups) * p * p) {
    return Status::failure(StatusCode::invalid_argument, "covariance size differs from n_groups * p * p");
  }
  for (std::int64_t k = 0; k < values.size; ++k) {
    if (!std::isfinite(values[k])) {
      return Status::failure(StatusCode::invalid_argument, "covariance input must be finite");
    }
  }
  std::vector<std::vector<std::int64_t>> rows_of_group(n_groups);
  for (std::int64_t i = 0; i < n; ++i) {
    const int g = group[i];
    if (g >= 0 && g < n_groups) rows_of_group[g].push_back(i);
  }

  auto body = [&](std::int64_t begin, std::int64_t end) {
    std::vector<double> column_mean(p);
    for (std::int64_t g = begin; g < end; ++g) {
      const std::vector<std::int64_t>& rows = rows_of_group[g];
      const std::int64_t n_obs = static_cast<std::int64_t>(rows.size());
      const std::int64_t offset = g * p * p;
      if (n_obs < 2) {
        for (std::int64_t k = 0; k < p * p; ++k) covariance[offset + k] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      // Column means with R's one-step refinement.
      for (std::int64_t j = 0; j < p; ++j) {
        const double* column = values.data + j * n;
        long double sum = 0.0L;
        for (std::int64_t row : rows) sum += column[row];
        long double mean = sum / static_cast<long double>(n_obs);
        if (std::isfinite(static_cast<double>(mean))) {
          sum = 0.0L;
          for (std::int64_t row : rows) sum += (column[row] - mean);
          mean = mean + sum / static_cast<long double>(n_obs);
        }
        column_mean[j] = static_cast<double>(mean);
      }
      // Lower triangle, mirrored.
      for (std::int64_t j = 0; j < p; ++j) {
        const double* column_j = values.data + j * n;
        for (std::int64_t l = 0; l <= j; ++l) {
          const double* column_l = values.data + l * n;
          long double sum = 0.0L;
          for (std::int64_t row : rows) {
            const double centred_j = column_j[row] - column_mean[j];
            const double centred_l = column_l[row] - column_mean[l];
            const double product = centred_j * centred_l;
            sum += product;
          }
          const double value = static_cast<double>(sum / static_cast<long double>(n_obs - 1));
          covariance[offset + l + j * p] = value;
          covariance[offset + j + l * p] = value;
        }
      }
    }
  };
  return parallel_for(n_groups, n_threads, 1, body, interrupted);
}

Status single_frame_statistics(const CscView& counts, Span<const double> n_count,
                               Span<const int> group, int n_groups, Span<double> focal_mean,
                               Span<double> within_ss, Span<double> global_mean, int n_threads,
                               const InterruptCheck& interrupted) {
  const std::int64_t n_cells = counts.n_rows;
  const std::int64_t n_genes = counts.n_cols;
  if (n_count.size != n_cells || group.size != n_cells) {
    return Status::failure(StatusCode::invalid_argument, "n_count and group must have one entry per cell");
  }
  if (focal_mean.size != static_cast<std::int64_t>(n_groups) * n_genes ||
      within_ss.size != static_cast<std::int64_t>(n_groups) * n_genes || global_mean.size != n_genes) {
    return Status::failure(StatusCode::invalid_argument, "output sizes do not match groups x genes");
  }
  for (std::int64_t i = 0; i < n_cells; ++i) {
    if (!std::isfinite(n_count[i]) || !(n_count[i] > 0.0)) {
      return Status::failure(StatusCode::invalid_argument,
                             "every cell needs a finite, positive total count (library size)");
    }
  }
  for (std::int64_t k = 0; k < counts.values.size; ++k) {
    if (!std::isfinite(counts.values[k])) {
      return Status::failure(StatusCode::invalid_argument, "counts must be finite (no NA, NaN or Inf)");
    }
  }
  std::vector<std::int64_t> group_size(n_groups, 0);
  for (std::int64_t i = 0; i < n_cells; ++i) {
    const int g = group[i];
    if (g >= 0 && g < n_groups) group_size[g] += 1;
  }

  auto body = [&](std::int64_t begin, std::int64_t end) {
    std::vector<long double> sum(n_groups);
    std::vector<std::int64_t> stored(n_groups);
    std::vector<double> mean(n_groups);
    std::vector<long double> squares(n_groups);
    for (std::int64_t gene = begin; gene < end; ++gene) {
      const int first = counts.column_pointer[gene];
      const int last = counts.column_pointer[gene + 1];
      std::fill(sum.begin(), sum.end(), 0.0L);
      std::fill(stored.begin(), stored.end(), 0);
      long double global_sum = 0.0L;
      // Pass 1: sums of y per group and overall.
      for (int k = first; k < last; ++k) {
        const int cell = counts.row_index[k];
        const double scale = 1e4 / n_count[cell];
        const double scaled = counts.values[k] * scale;
        const double y = std::log1p(scaled);
        global_sum += y;
        const int g = group[cell];
        if (g < 0 || g >= n_groups) continue;
        sum[g] += y;
        stored[g] += 1;
      }
      global_mean[gene] = static_cast<double>(global_sum / static_cast<long double>(n_cells));
      for (int g = 0; g < n_groups; ++g) {
        mean[g] = static_cast<double>(sum[g] / static_cast<long double>(group_size[g]));
      }
      // Pass 2: centred squares over stored entries, plus the zeros' contribution.
      std::fill(squares.begin(), squares.end(), 0.0L);
      for (int k = first; k < last; ++k) {
        const int cell = counts.row_index[k];
        const int g = group[cell];
        if (g < 0 || g >= n_groups) continue;
        const double scale = 1e4 / n_count[cell];
        const double scaled = counts.values[k] * scale;
        const double y = std::log1p(scaled);
        const double deviation = y - mean[g];
        const double deviation_sq = deviation * deviation;
        squares[g] += deviation_sq;
      }
      for (int g = 0; g < n_groups; ++g) {
        const std::int64_t slot = g + gene * static_cast<std::int64_t>(n_groups);
        focal_mean[slot] = mean[g];
        if (group_size[g] == 0) {
          within_ss[slot] = 0.0;
          continue;
        }
        const double zero_deviation_sq = mean[g] * mean[g];
        const long double zeros = static_cast<long double>(group_size[g] - stored[g]);
        within_ss[slot] = static_cast<double>(squares[g] + zeros * zero_deviation_sq);
      }
    }
  };
  return parallel_for(n_genes, n_threads, 16, body, interrupted);
}

}  // namespace pace
