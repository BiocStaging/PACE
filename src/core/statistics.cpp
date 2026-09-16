// statistics.cpp -- implementation of statistics.hpp.
#include "statistics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "r_summaries.hpp"
#include "thread_pool.hpp"

namespace pace {

bool all_finite(Span<const double> values) {
  for (std::int64_t k = 0; k < values.size; ++k) {
    if (!std::isfinite(values[k])) return false;
  }
  return true;
}

Status dense_group_moments(Span<const double> values, std::int64_t n, std::int64_t p,
                           Span<const int> group, int n_groups, Span<double> mean,
                           Span<double> variance, Span<int> n_used, bool* any_nonzero,
                           int n_threads, const InterruptCheck& interrupted) {
  const bool all_zero = values.size == 0;
  if (!all_zero && values.size != n * p) {
    return Status::failure(StatusCode::invalid_argument, "values must be empty or n * p");
  }
  if (group.size != n) {
    return Status::failure(StatusCode::invalid_argument, "group must have one entry per row");
  }
  const std::int64_t slots = static_cast<std::int64_t>(n_groups) * p;
  if ((mean.size != 0 && mean.size != slots) || (variance.size != 0 && variance.size != slots) ||
      (n_used.size != 0 && n_used.size != slots)) {
    return Status::failure(StatusCode::invalid_argument, "output sizes must be empty or n_groups * p");
  }

  std::vector<std::vector<std::int64_t>> rows_of_group(n_groups);
  for (std::int64_t i = 0; i < n; ++i) {
    const int g = group[i];
    if (g >= 0 && g < n_groups) rows_of_group[g].push_back(i);
  }
  // One flag per column, so the reduction does not depend on the thread count.
  std::vector<char> column_nonzero(static_cast<std::size_t>(p), 0);

  auto body = [&](std::int64_t begin, std::int64_t end) {
    std::vector<double> zeros;
    for (std::int64_t j = begin; j < end; ++j) {
      const double* column = nullptr;
      if (all_zero) {
        if (zeros.empty()) zeros.assign(static_cast<std::size_t>(n), 0.0);
        column = zeros.data();
      } else {
        column = values.data + j * n;
      }
      if (any_nonzero != nullptr && !all_zero) {
        for (std::int64_t i = 0; i < n; ++i) {
          const double value = column[i];
          if (!std::isnan(value) && value != 0.0) {
            column_nonzero[static_cast<std::size_t>(j)] = 1;
            break;
          }
        }
      }
      for (int g = 0; g < n_groups; ++g) {
        const std::int64_t slot = g + j * static_cast<std::int64_t>(n_groups);
        std::int64_t used = 0;
        const double group_mean = r_mean(column, rows_of_group[g], true, &used);
        if (mean.size != 0) mean[slot] = group_mean;
        if (variance.size != 0) {
          variance[slot] = r_variance(column, rows_of_group[g], group_mean, used, true);
        }
        if (n_used.size != 0) n_used[slot] = static_cast<int>(used);
      }
    }
  };
  const Status status = parallel_for(p, n_threads, 8, body, interrupted);
  if (!status.is_ok()) return status;
  if (any_nonzero != nullptr) {
    *any_nonzero = false;
    for (std::int64_t j = 0; j < p; ++j) {
      if (column_nonzero[static_cast<std::size_t>(j)] != 0) *any_nonzero = true;
    }
  }
  return Status::success();
}

Status final_pass_statistics(Span<const double> eta, std::int64_t n, std::int64_t n_genes,
                             Span<const double> offset, Span<const double> ambient,
                             Span<const double> rho, Span<const int> group, int n_groups,
                             Span<double> mu_group_sum, Span<double> toff_variance,
                             Span<double> mu_column_sum, Span<double> spill_row_sum,
                             Span<double> total_row_sum, Span<double> mu_out,
                             Span<double> toff_out, bool* any_nonzero,
                             const InterruptCheck& interrupted) {
  const std::int64_t cells_times_genes = n * n_genes;
  if (eta.size != cells_times_genes || ambient.size != cells_times_genes) {
    return Status::failure(StatusCode::invalid_argument, "eta and ambient must be n * n_genes");
  }
  if (offset.size != n || rho.size != n || group.size != n) {
    return Status::failure(StatusCode::invalid_argument, "offset, rho and group must have one entry per cell");
  }
  const std::int64_t slots = static_cast<std::int64_t>(n_groups) * n_genes;
  if (mu_group_sum.size != slots || toff_variance.size != slots || mu_column_sum.size != n_genes) {
    return Status::failure(StatusCode::invalid_argument, "per-group outputs must be n_groups * n_genes");
  }
  if (spill_row_sum.size != n || total_row_sum.size != n) {
    return Status::failure(StatusCode::invalid_argument, "row-sum outputs must have one entry per cell");
  }
  if ((mu_out.size != 0 && mu_out.size != cells_times_genes) ||
      (toff_out.size != 0 && toff_out.size != cells_times_genes)) {
    return Status::failure(StatusCode::invalid_argument, "mu and toff outputs must be empty or n * n_genes");
  }
  for (std::int64_t i = 0; i < n; ++i) {
    if (!std::isfinite(offset[i]) || !std::isfinite(rho[i])) {
      return Status::failure(StatusCode::invalid_argument, "offset and rho must be finite");
    }
  }
  for (std::int64_t k = 0; k < cells_times_genes; ++k) {
    if (!std::isfinite(eta[k]) || !std::isfinite(ambient[k])) {
      return Status::failure(StatusCode::invalid_argument, "eta and the ambient field must be finite");
    }
  }

  std::vector<std::vector<std::int64_t>> rows_of_group(n_groups);
  for (std::int64_t i = 0; i < n; ++i) {
    const int g = group[i];
    if (g >= 0 && g < n_groups) rows_of_group[g].push_back(i);
  }
  std::vector<long double> spill_accumulator(static_cast<std::size_t>(n), 0.0L);
  std::vector<long double> total_accumulator(static_cast<std::size_t>(n), 0.0L);
  std::vector<double> mu_column(static_cast<std::size_t>(n));
  std::vector<double> toff_column(static_cast<std::size_t>(n));
  if (any_nonzero != nullptr) *any_nonzero = false;

  for (std::int64_t j = 0; j < n_genes; ++j) {
    if ((j % 16) == 0 && interrupted && interrupted()) {
      return Status::failure(StatusCode::interrupted, "interrupted");
    }
    const double* eta_column = eta.data + j * n;
    const double* ambient_column = ambient.data + j * n;
    long double column_sum = 0.0L;
    for (std::int64_t i = 0; i < n; ++i) {
      double mu_bio = std::exp(eta_column[i] + offset[i]);
      if (!(mu_bio > 1e-6)) mu_bio = 1e-6;
      double spill = ambient_column[i] * rho[i];
      if (!(spill > 0.0)) spill = 0.0;
      double mu = mu_bio + spill;
      if (!(mu > 1e-6)) mu = 1e-6;
      const double toff = std::log1p(spill / mu_bio);
      mu_column[static_cast<std::size_t>(i)] = mu;
      toff_column[static_cast<std::size_t>(i)] = toff;
      column_sum += mu;
      spill_accumulator[static_cast<std::size_t>(i)] += spill;
      total_accumulator[static_cast<std::size_t>(i)] += mu;
      if (any_nonzero != nullptr && toff != 0.0 && !std::isnan(toff)) *any_nonzero = true;
    }
    mu_column_sum[j] = static_cast<double>(column_sum);
    for (int g = 0; g < n_groups; ++g) {
      const std::vector<std::int64_t>& rows = rows_of_group[g];
      const std::int64_t slot = g + j * static_cast<std::int64_t>(n_groups);
      long double group_sum = 0.0L;
      for (std::int64_t row : rows) group_sum += mu_column[static_cast<std::size_t>(row)];
      mu_group_sum[slot] = static_cast<double>(group_sum);
      std::int64_t used = 0;
      const double group_mean = r_mean(toff_column.data(), rows, false, &used);
      toff_variance[slot] = r_variance(toff_column.data(), rows, group_mean, used, false);
    }
    if (mu_out.size != 0) {
      std::copy(mu_column.begin(), mu_column.end(), mu_out.data + j * n);
    }
    if (toff_out.size != 0) {
      std::copy(toff_column.begin(), toff_column.end(), toff_out.data + j * n);
    }
  }
  for (std::int64_t i = 0; i < n; ++i) {
    spill_row_sum[i] = static_cast<double>(spill_accumulator[static_cast<std::size_t>(i)]);
    total_row_sum[i] = static_cast<double>(total_accumulator[static_cast<std::size_t>(i)]);
  }
  return Status::success();
}

}  // namespace pace
