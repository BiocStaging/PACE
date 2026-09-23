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

// One gene column of a CSC matrix, expanded over the rows [lo, hi) only, into
// `out[0 .. hi - lo)`. A worker that owns a slice of the cells must not pay the
// O(n) assign() that expand_column() does, once per gene, per worker. Row
// indices within a dgCMatrix column are sorted, so the slice is one binary
// search away.
void expand_column_range(const GeneBlock& block, std::int64_t gene, std::int64_t lo,
                         std::int64_t hi, std::vector<double>& out) {
  out.assign(static_cast<std::size_t>(hi - lo), 0.0);
  if (block.matrix.column_pointer.size == 0) return;
  const std::int64_t column = block.first_gene + gene;
  const int first = block.matrix.column_pointer[column];
  const int last = block.matrix.column_pointer[column + 1];
  const int* rows = block.matrix.row_index.data;
  const int* begin = std::lower_bound(rows + first, rows + last, static_cast<int>(lo));
  const int* end = std::lower_bound(begin, rows + last, static_cast<int>(hi));
  for (const int* it = begin; it != end; ++it) {
    const std::int64_t k = it - rows;
    out[static_cast<std::size_t>(*it - lo)] = block.matrix.values[k];
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

// The weight one (cell, gene) puts on the ambient field: the NB precision at the
// total mean, times the ambient value. Shared by the serial and the parallel
// paths below so that the arithmetic cannot drift between them.
inline double rho_weighted_ambient(bool seed_iteration, double count, double eta_value,
                                   double offset_value, double ambient_value, double rho_value,
                                   bool nb2, double alpha_gene) {
  double mu_total;
  if (seed_iteration) {
    // The solver's seed: mu = max(y, 0.5), with no contamination part yet.
    mu_total = floor_at(count, 0.5);
  } else {
    const double mu_bio = floor_at(std::exp(eta_value + offset_value), 1e-6);
    const double spill = floor_at(ambient_value * rho_value, 0.0);
    mu_total = floor_at(mu_bio + spill, 1e-8);
  }
  const double precision = nb2 ? 1 / (mu_total * (1 + mu_total * alpha_gene))
                               : (1 / mu_total) / (1 + alpha_gene);
  return precision * ambient_value;
}


}  // namespace

Status working_response(Span<const double> eta, const GeneBlock& counts, const GeneBlock& ambient,
                        Span<const double> offset, Span<const double> rho,
                        Span<const double> alpha, Span<const double> sample_weight, bool nb2,
                        bool gaussian, bool seed_iteration, std::int64_t n, std::int64_t n_genes,
                        Span<double> z, Span<double> w, Span<double> colsum_w, int n_threads,
                        const InterruptCheck& interrupted) {
  const std::int64_t cells_times_genes = n * n_genes;
  // The Gaussian path reads neither eta nor the ambient block, so it is exempt
  // from the checks on them; `alpha` carries sigma2_g and must still be sized.
  if (!gaussian && !seed_iteration && eta.size != cells_times_genes) {
    return Status::failure(StatusCode::invalid_argument, "eta must be n * n_genes");
  }
  if (gaussian) {
    if (alpha.size != n_genes) {
      return Status::failure(StatusCode::invalid_argument, "sigma2 must have one entry per gene");
    }
  } else if (offset.size != n || rho.size != n || alpha.size != n_genes) {
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
  if (!gaussian && !seed_iteration && !block_covers(ambient, n, n_genes)) {
    return Status::failure(StatusCode::invalid_argument, "the ambient block does not cover the genes");
  }
  if (!gaussian) {
    for (std::int64_t i = 0; i < n; ++i) {
      if (!std::isfinite(offset[i]) || !std::isfinite(rho[i])) {
        return Status::failure(StatusCode::invalid_argument, "offset and rho must be finite");
      }
    }
  } else {
    for (std::int64_t j = 0; j < n_genes; ++j) {
      if (!(alpha[j] > 0) || !std::isfinite(alpha[j])) {
        return Status::failure(StatusCode::invalid_argument, "sigma2 must be finite and positive");
      }
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
      if (gaussian) {
        // Identity link: z is the raw intensity and w is the per-gene scalar
        // 1 / sigma2_g, so neither depends on eta. The R engine builds the same
        // weight with rep(1 / sigma2[genes], each = n); a scalar written n times
        // and a scalar summed n times agree with it to the bit, and colsum_w
        // accumulates in long double exactly as R's colSums does.
        const double inverse_sigma2 = 1.0 / alpha[j];
        long double weight_sum = 0.0L;
        for (std::int64_t i = 0; i < n; ++i) {
          double weight = inverse_sigma2;
          if (sample_weight.size != 0) weight = weight * sample_weight[i];
          z[i + j * n] = y_column[static_cast<std::size_t>(i)];
          w[i + j * n] = weight;
          if (!std::isnan(weight)) weight_sum += weight;
        }
        colsum_w[j] = static_cast<double>(weight_sum);
        continue;
      }
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
                      int n_threads, const InterruptCheck& interrupted) {
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

  // ---- the convergence metric, serially, in gene-then-cell order ----
  // Kept serial so delta_sum, a long double summed across both axes, associates
  // exactly as it always has. It carries no exp() and no sparse expansion unless
  // the previous eta has to be seeded from the counts, so it is the cheap half.
  if (want_delta) {
    for (std::int64_t j = 0; j < n_genes; ++j) {
      if ((j % 8) == 0 && interrupted && interrupted()) {
        return Status::failure(StatusCode::interrupted, "interrupted");
      }
      if (!have_previous) expand_column(counts, j, n, y_column);
      const double* eta_column = eta.data + j * n;
      const double* prev_column = have_previous ? prev_eta.data + j * n : nullptr;
      for (std::int64_t i = 0; i < n; ++i) {
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
    }
  }

  // Without a convergence metric there is nothing accumulated across cells, only
  // into them: num[i] and den[i] are per-cell, summed over genes. So the CELLS
  // can be split across workers while the GENES stay in order, and every cell is
  // touched by exactly one worker. Each cell's long double sum then runs over the
  // same genes in the same order it did serially, which is what makes this
  // bit-identical rather than merely deterministic.
  //
  // The parallelism sits ABOVE the gene loop on purpose: parallel_for starts and
  // joins fresh threads on every call, so dispatching once per gene would spend
  // more on thread creation than the loop costs (about 300 us a call, against
  // 14,170 calls for a 1,090-gene cohort over 13 iterations).
  //
  // The convergence metric is taken FIRST, in its own serial pass just above,
  // because its delta_sum is a long double accumulated across both cells and
  // genes: splitting the cells would reassociate it. Separating the two costs a
  // second traversal but keeps every number identical, and the cheap half is the
  // one left serial -- the delta is a subtraction, an abs and a divide, while the
  // half that moves to the workers carries the exp() in mu_bio.
  if (n_threads > 1) {
    const std::int64_t cells_per_block = 8192;   // fixed: the partition must not
                                                 // depend on the worker count
    const Status status = parallel_for(
        n, n_threads, cells_per_block,
        [&](std::int64_t lo, std::int64_t hi) {
          std::vector<double> y_slice;
          std::vector<double> ambient_slice;
          for (std::int64_t j = 0; j < n_genes; ++j) {
            expand_column_range(counts, j, lo, hi, y_slice);
            expand_column_range(ambient, j, lo, hi, ambient_slice);
            const double* eta_column = seed_iteration ? nullptr : eta.data + j * n;
            const double alpha_gene = alpha[j];
            // counts.first_gene, NOT ambient.first_gene. Both are the global gene
            // index when the ambient field is a slice of a cached n x G product,
            // which is why using either worked. They are not the same thing when
            // the ambient block is built per chunk: that block is chunk-local and
            // starts at zero, so deriving the mask column from it reads the wrong
            // anchors for every chunk after the first. The counts block is always
            // a slice of the full panel, so its offset is the global index.
            const std::int64_t mask_column = counts.first_gene + j;
            for (std::int64_t i = lo; i < hi; ++i) {
              const std::size_t slot = static_cast<std::size_t>(i - lo);
              const double ambient_value = ambient_slice[slot];
              double weighted_ambient = rho_weighted_ambient(
                  seed_iteration, y_slice[slot], seed_iteration ? 0.0 : eta_column[i], offset[i],
                  ambient_value, rho[i], nb2, alpha_gene);
              if (mask.size != 0) {
                const std::int64_t row = mask_index.size != 0 ? (mask_index[i] - 1) : i;
                weighted_ambient = weighted_ambient * mask[row + mask_column * n_mask_rows];
              }
              const double numerator = weighted_ambient * y_slice[slot];
              const double denominator = weighted_ambient * ambient_value;
              if (!std::isnan(numerator)) num_accumulator[static_cast<std::size_t>(i)] += numerator;
              if (!std::isnan(denominator)) den_accumulator[static_cast<std::size_t>(i)] += denominator;
            }
          }
        },
        interrupted);
    if (!status.is_ok()) return status;
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

  for (std::int64_t j = 0; j < n_genes; ++j) {
    if ((j % 8) == 0 && interrupted && interrupted()) {
      return Status::failure(StatusCode::interrupted, "interrupted");
    }
    expand_column(counts, j, n, y_column);
    expand_column(ambient, j, n, ambient_column);
    const double* eta_column = seed_iteration ? nullptr : eta.data + j * n;
    const double alpha_gene = alpha[j];
    const std::int64_t mask_column = counts.first_gene + j;   // see the note in the parallel path
    for (std::int64_t i = 0; i < n; ++i) {
      const double ambient_value = ambient_column[static_cast<std::size_t>(i)];
      double weighted_ambient = rho_weighted_ambient(
          seed_iteration, y_column[static_cast<std::size_t>(i)],
          seed_iteration ? 0.0 : eta_column[i], offset[i], ambient_value, rho[i], nb2, alpha_gene);
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
