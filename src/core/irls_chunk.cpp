// irls_chunk.cpp -- implementation of irls_chunk.hpp.
#include "fp_no_contract.hpp"  // must precede the arithmetic below

#include "irls_chunk.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "r_summaries.hpp"
#include "thread_pool.hpp"

namespace pace {
namespace {

// One gene column of a CSC matrix, expanded into a dense scratch vector.
void expand_column(const GeneBlock& block, std::int64_t gene, std::int64_t n,
                   std::vector<double>& out) {
  out.assign(static_cast<std::size_t>(n), 0.0);
  if (block.matrix.column_pointer.size == 0) return;
  const std::int64_t column = block.first_gene + gene;
  for (int k = block.matrix.column_pointer[column]; k < block.matrix.column_pointer[column + 1];
       ++k) {
    out[static_cast<std::size_t>(block.matrix.row_index[k])] = block.matrix.values[k];
  }
}

bool block_covers(const GeneBlock& block, std::int64_t n, std::int64_t n_genes) {
  if (block.matrix.column_pointer.size == 0) return false;
  return block.matrix.n_rows == n && block.first_gene + n_genes <= block.matrix.n_cols;
}

// R's pmax(value, floor) for one value, then the value itself: a missing value
// stays missing, as pmax() leaves NA.
double floor_at(double value, double floor_value) {
  if (std::isnan(value)) return value;
  return value > floor_value ? value : floor_value;
}

}  // namespace

Status working_response(Span<const double> eta, const GeneBlock& counts, const GeneBlock& ambient,
                        Span<const double> offset, Span<const double> rho,
                        Span<const double> alpha, Span<const double> sample_weight, bool nb2,
                        bool seed_iteration, std::int64_t n, std::int64_t n_genes, Span<double> z,
                        Span<double> w, Span<double> colsum_w, int n_threads,
                        const InterruptCheck& interrupted) {
  const std::int64_t cells_times_genes = n * n_genes;
  if (!seed_iteration && eta.size != cells_times_genes) {
    return Status::failure(StatusCode::invalid_argument, "eta must be n * n_genes");
  }
  if (offset.size != n || rho.size != n || alpha.size != n_genes) {
    return Status::failure(StatusCode::invalid_argument, "offset, rho or alpha have the wrong size");
  }
  if (sample_weight.size != 0 && sample_weight.size != n) {
    return Status::failure(StatusCode::invalid_argument, "sample_weight must be empty or one per cell");
  }
  if (z.size != cells_times_genes || w.size != cells_times_genes || colsum_w.size != n_genes) {
    return Status::failure(StatusCode::invalid_argument, "outputs must be n * n_genes");
  }
  if (!block_covers(counts, n, n_genes)) {
    return Status::failure(StatusCode::invalid_argument, "the count block does not cover the genes");
  }
  if (!seed_iteration && !block_covers(ambient, n, n_genes)) {
    return Status::failure(StatusCode::invalid_argument, "the ambient block does not cover the genes");
  }
  for (std::int64_t i = 0; i < n; ++i) {
    if (!std::isfinite(offset[i]) || !std::isfinite(rho[i])) {
      return Status::failure(StatusCode::invalid_argument, "offset and rho must be finite");
    }
  }

  // The scratch is thread_local, not per call. With one gene per block (below)
  // `body` runs once per gene, and a fresh pair of n-length vectors per gene
  // would be 19.6 MB of malloc and first-touch faults each at 1.2M cells. Held
  // per worker instead, expand_column's assign() reuses the capacity and the
  // allocation happens once per thread for the life of the pool.
  auto body = [&](std::int64_t begin, std::int64_t end) {
    static thread_local std::vector<double> y_column;
    static thread_local std::vector<double> ambient_column;
    for (std::int64_t j = begin; j < end; ++j) {
      expand_column(counts, j, n, y_column);
      if (!seed_iteration) expand_column(ambient, j, n, ambient_column);
      const double* eta_column = seed_iteration ? nullptr : eta.data + j * n;
      const double alpha_gene = alpha[j];
      long double weight_sum = 0.0L;
      for (std::int64_t i = 0; i < n; ++i) {
        const double y = y_column[static_cast<std::size_t>(i)];
        double mu_bio;
        double mu;
        if (seed_iteration) {
          mu_bio = floor_at(y, 0.5);
          mu = mu_bio;
        } else {
          mu_bio = floor_at(std::exp(eta_column[i] + offset[i]), 1e-6);
          const double spill = floor_at(ambient_column[static_cast<std::size_t>(i)] * rho[i], 0.0);
          mu = floor_at(mu_bio + spill, 1e-6);
        }
        const double working_eta = std::log(mu_bio) - offset[i];
        double weight;
        if (nb2) {
          weight = (mu_bio * mu_bio) / (mu * (1 + mu * alpha_gene));
        } else {
          weight = ((mu_bio * mu_bio) / mu) / (1 + alpha_gene);
        }
        if (sample_weight.size != 0) weight = weight * sample_weight[i];
        z[i + j * n] = working_eta + (y - mu) / mu_bio;
        w[i + j * n] = weight;
        if (!std::isnan(weight)) weight_sum += weight;
      }
      colsum_w[j] = static_cast<double>(weight_sum);
    }
  };
  // One gene per block, not four. Four blocks meant worker_cap = min(threads, 4)
  // in thread_pool, so this pass used at most FOUR threads however many were
  // asked for -- the caller hands it sub-blocks of 16 genes (bindings.cpp), so
  // four genes a block left twelve of those sixteen unable to find a worker.
  // Per gene, the ceiling becomes the sub-block width instead, and the
  // partition still does not depend on the worker count, so neither do results:
  // gene j writes only z[, j], w[, j] and colsum_w[j], and the long double
  // weight_sum never crosses a gene boundary.
  return parallel_for(n_genes, n_threads, 1, body, interrupted);
}

Status rho_accumulate(Span<const double> eta, Span<const double> prev_eta, const GeneBlock& counts,
                      const GeneBlock& ambient, Span<const double> offset,
                      Span<const double> rho, Span<const double> alpha, Span<const double> mask,
                      Span<const int> mask_index, std::int64_t n_mask_rows, bool nb2,
                      bool seed_iteration, bool seed_previous, std::int64_t n,
                      std::int64_t n_genes, Span<double> num,
                      Span<double> den, double* rel_delta_max, double* rel_delta_sum,
                      std::int64_t* n_finite, std::int64_t* n_nonfinite, Span<double> tail_counts,
                      const InterruptCheck& interrupted) {
  const std::int64_t cells_times_genes = n * n_genes;
  const bool have_previous = prev_eta.size == cells_times_genes;
  const bool want_delta = have_previous || seed_previous;
  if ((!seed_iteration && eta.size != cells_times_genes) ||
      (prev_eta.size != 0 && !have_previous)) {
    return Status::failure(StatusCode::invalid_argument, "eta and prev_eta must be n * n_genes");
  }
  if (offset.size != n || rho.size != n || alpha.size != n_genes) {
    return Status::failure(StatusCode::invalid_argument, "offset, rho or alpha have the wrong size");
  }
  if (num.size != n || den.size != n) {
    return Status::failure(StatusCode::invalid_argument, "num and den must have one entry per cell");
  }
  if (!block_covers(counts, n, n_genes) || !block_covers(ambient, n, n_genes)) {
    return Status::failure(StatusCode::invalid_argument, "a block does not cover the genes");
  }
  if (mask.size != 0 && mask_index.size != 0 && mask_index.size != n) {
    return Status::failure(StatusCode::invalid_argument, "the mask index must have one entry per cell");
  }

  // Row sums in long double across the whole chunk, in the column order R used.
  std::vector<long double> num_accumulator(static_cast<std::size_t>(n), 0.0L);
  std::vector<long double> den_accumulator(static_cast<std::size_t>(n), 0.0L);
  std::vector<double> y_column;
  std::vector<double> ambient_column;
  long double delta_sum = 0.0L;
  double delta_max = *rel_delta_max;
  std::int64_t finite_count = 0;
  std::int64_t nonfinite_count = 0;

  for (std::int64_t j = 0; j < n_genes; ++j) {
    if ((j % 8) == 0 && interrupted && interrupted()) {
      return Status::failure(StatusCode::interrupted, "interrupted");
    }
    expand_column(counts, j, n, y_column);
    expand_column(ambient, j, n, ambient_column);
    const double* eta_column = seed_iteration ? nullptr : eta.data + j * n;
    const double* prev_column = have_previous ? prev_eta.data + j * n : nullptr;
    const double alpha_gene = alpha[j];
    const std::int64_t mask_column = ambient.first_gene + j;
    for (std::int64_t i = 0; i < n; ++i) {
      if (want_delta) {
        const double previous =
            have_previous ? prev_column[i]
                          : std::log(floor_at(y_column[static_cast<std::size_t>(i)], 0.5)) - offset[i];
        const double delta = std::abs(eta_column[i] - previous) / floor_at(std::abs(previous), 1e-3);
        if (std::isfinite(delta)) {
          if (!(delta <= delta_max)) delta_max = delta;
          delta_sum += delta;
          finite_count += 1;
          if (tail_counts.size == 4) {
            if (delta > 0.01) tail_counts[0] += 1;
            if (delta > 0.05) tail_counts[1] += 1;
            if (delta > 0.1) tail_counts[2] += 1;
            if (delta > 1.0) tail_counts[3] += 1;
          }
        } else {
          // max(x, na.rm = TRUE) skips NaN but keeps an infinite delta.
          if (!std::isnan(delta) && !(delta <= delta_max)) delta_max = delta;
          nonfinite_count += 1;
        }
      }
      const double ambient_value = ambient_column[static_cast<std::size_t>(i)];
      double mu_total;
      if (seed_iteration) {
        // The solver's seed: mu = max(y, 0.5), with no contamination part yet.
        mu_total = floor_at(y_column[static_cast<std::size_t>(i)], 0.5);
      } else {
        const double mu_bio = floor_at(std::exp(eta_column[i] + offset[i]), 1e-6);
        const double spill = floor_at(ambient_value * rho[i], 0.0);
        mu_total = floor_at(mu_bio + spill, 1e-8);
      }
      double precision;
      if (nb2) {
        precision = 1 / (mu_total * (1 + mu_total * alpha_gene));
      } else {
        precision = (1 / mu_total) / (1 + alpha_gene);
      }
      double weighted_ambient = precision * ambient_value;
      if (mask.size != 0) {
        const std::int64_t row = mask_index.size != 0 ? (mask_index[i] - 1) : i;
        weighted_ambient = weighted_ambient * mask[row + mask_column * n_mask_rows];
      }
      const double numerator = weighted_ambient * y_column[static_cast<std::size_t>(i)];
      const double denominator = weighted_ambient * ambient_value;
      if (!std::isnan(numerator)) num_accumulator[static_cast<std::size_t>(i)] += numerator;
      if (!std::isnan(denominator)) den_accumulator[static_cast<std::size_t>(i)] += denominator;
    }
  }
  for (std::int64_t i = 0; i < n; ++i) {
    num[i] += static_cast<double>(num_accumulator[static_cast<std::size_t>(i)]);
    den[i] += static_cast<double>(den_accumulator[static_cast<std::size_t>(i)]);
  }
  *rel_delta_max = delta_max;
  *rel_delta_sum += static_cast<double>(delta_sum);
  *n_finite += finite_count;
  *n_nonfinite += nonfinite_count;
  return Status::success();
}

Status fitted_mean_column(Span<const double> eta_column, const GeneBlock& counts,
                          const GeneBlock& ambient, Span<const double> offset,
                          Span<const double> rho, std::int64_t n, std::int64_t gene,
                          Span<double> mu, Span<double> y) {
  if (eta_column.size != n || offset.size != n || rho.size != n || mu.size != n || y.size != n) {
    return Status::failure(StatusCode::invalid_argument, "every vector must have n entries");
  }
  std::vector<double> y_column;
  std::vector<double> ambient_column;
  expand_column(counts, gene, n, y_column);
  expand_column(ambient, gene, n, ambient_column);
  for (std::int64_t i = 0; i < n; ++i) {
    const double mu_bio = floor_at(std::exp(eta_column[i] + offset[i]), 1e-6);
    const double spill = floor_at(ambient_column[static_cast<std::size_t>(i)] * rho[i], 0.0);
    mu[i] = floor_at(mu_bio + spill, 1e-6);
    y[i] = y_column[static_cast<std::size_t>(i)];
  }
  return Status::success();
}

Status rho_shrink(Span<const double> num, Span<const double> den, std::int64_t n, Span<double> rho,
                  std::int64_t* n_nonfinite, double* rho_bar, double* den_quantile) {
  if (num.size != n || den.size != n || rho.size != n) {
    return Status::failure(StatusCode::invalid_argument, "every vector must have n entries");
  }
  std::vector<double> raw(static_cast<std::size_t>(n), 0.0);
  std::vector<double> weight(static_cast<std::size_t>(n), 0.0);
  std::int64_t bad = 0;
  for (std::int64_t i = 0; i < n; ++i) {
    const bool finite_cell = std::isfinite(den[i]) && std::isfinite(num[i]);
    if (!finite_cell) bad += 1;
    double ratio = 0.0;
    if (finite_cell && den[i] > 1e-12) ratio = num[i] / den[i];
    raw[static_cast<std::size_t>(i)] = ratio > 0.0 ? ratio : 0.0;
    if (finite_cell && den[i] > 0.0) weight[static_cast<std::size_t>(i)] = den[i];
  }
  long double weighted = 0.0L;
  long double weight_total = 0.0L;
  for (std::int64_t i = 0; i < n; ++i) {
    // the product is rounded to double before it enters the sum, as R's
    // sum(den_w * rho_raw) does
    const double product = weight[static_cast<std::size_t>(i)] * raw[static_cast<std::size_t>(i)];
    weighted += product;
    weight_total += weight[static_cast<std::size_t>(i)];
  }
  const double weight_sum = static_cast<double>(weight_total);
  const double bar = static_cast<double>(weighted) / (weight_sum > 1e-12 ? weight_sum : 1e-12);
  *rho_bar = bar;
  *n_nonfinite = bad;

  std::vector<double> positive;
  positive.reserve(static_cast<std::size_t>(n));
  for (std::int64_t i = 0; i < n; ++i) {
    if (weight[static_cast<std::size_t>(i)] > 1e-12) positive.push_back(weight[static_cast<std::size_t>(i)]);
  }
  if (positive.empty()) {
    for (std::int64_t i = 0; i < n; ++i) rho[i] = 0.0;
    *den_quantile = std::numeric_limits<double>::quiet_NaN();
    return Status::success();
  }
  // R's quantile(x, 0.10) with the default type 7.
  std::sort(positive.begin(), positive.end());
  const double index = 1 + (static_cast<double>(positive.size()) - 1) * 0.10;
  const double low = std::floor(index);
  const double high = std::ceil(index);
  const double at_low = positive[static_cast<std::size_t>(low) - 1];
  const double at_high = positive[static_cast<std::size_t>(high) - 1];
  double den0 = at_low;
  if (index > low && at_high != at_low) {
    const double h = index - low;
    den0 = (1 - h) * at_low + h * at_high;
  }
  *den_quantile = den0;
  for (std::int64_t i = 0; i < n; ++i) {
    const double value = (weight[static_cast<std::size_t>(i)] * raw[static_cast<std::size_t>(i)] +
                          den0 * bar) /
                         (weight[static_cast<std::size_t>(i)] + den0);
    rho[i] = floor_at(value, 0.0);
  }
  return Status::success();
}

}  // namespace pace
