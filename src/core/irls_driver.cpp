// irls_driver.cpp -- implementation of irls_driver.hpp.
#include "fp_no_contract.hpp"  // must precede the arithmetic below

#include "irls_driver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <limits>
#include <vector>

#include "hyperparameters.hpp"
#include "irls_chunk.hpp"
#include "linear_predictor.hpp"
#include "r_summaries.hpp"
#include "statistics.hpp"

namespace pace {
namespace {

const double infinity = std::numeric_limits<double>::infinity();
const double not_a_number = std::numeric_limits<double>::quiet_NaN();

Span<const double> empty_doubles() { return Span<const double>(nullptr, 0); }

Span<double> empty_output() { return Span<double>(nullptr, 0); }

// One gene column range of a CSC block, as the pass loops hand it on.
GeneBlock shifted_block(const GeneBlock& block, std::int64_t first_gene) {
  GeneBlock out = block;
  out.first_gene = block.first_gene + first_gene;
  return out;
}

// sprintf into a std::string, so the progress lines are written the way the R
// `cat(sprintf(...))` they replace wrote them.
std::string format_text(const char* format, ...) {
  char buffer[1024];
  std::va_list arguments;
  va_start(arguments, format);
  const int written = std::vsnprintf(buffer, sizeof(buffer), format, arguments);
  va_end(arguments);
  if (written < 0) return std::string();
  return std::string(buffer, static_cast<std::size_t>(
                                 written < static_cast<int>(sizeof(buffer))
                                     ? written
                                     : static_cast<int>(sizeof(buffer)) - 1));
}

// R's median(x) with na.rm = FALSE: a single missing value makes the whole
// median missing, which is why r_median() (na.rm = TRUE) cannot be used for the
// progress lines.
double r_median_strict(const double* values, std::int64_t size) {
  if (size == 0) return not_a_number;
  std::vector<double> scratch(static_cast<std::size_t>(size));
  for (std::int64_t i = 0; i < size; ++i) {
    if (std::isnan(values[i])) return not_a_number;
    scratch[static_cast<std::size_t>(i)] = values[i];
  }
  return r_median(scratch);
}

// R's quantile(x, probability) with the default type 7.
double r_quantile_type7(const double* values, std::int64_t size, double probability) {
  if (size == 0) return not_a_number;
  std::vector<double> scratch(static_cast<std::size_t>(size));
  for (std::int64_t i = 0; i < size; ++i) {
    if (std::isnan(values[i])) return not_a_number;
    scratch[static_cast<std::size_t>(i)] = values[i];
  }
  std::sort(scratch.begin(), scratch.end());
  const double index = 1 + (static_cast<double>(size) - 1) * probability;
  const double low = std::floor(index);
  const double high = std::ceil(index);
  const double at_low = scratch[static_cast<std::size_t>(low) - 1];
  const double at_high = scratch[static_cast<std::size_t>(high) - 1];
  if (index > low && at_high != at_low) {
    const double h = index - low;
    return (1 - h) * at_low + h * at_high;
  }
  return at_low;
}

// R's min()/max() over a double vector: a missing value propagates.
double r_min(const double* values, std::int64_t size) {
  double smallest = infinity;
  for (std::int64_t i = 0; i < size; ++i) {
    if (std::isnan(values[i])) return not_a_number;
    if (values[i] < smallest) smallest = values[i];
  }
  return smallest;
}

double r_max(const double* values, std::int64_t size) {
  double largest = -infinity;
  for (std::int64_t i = 0; i < size; ++i) {
    if (std::isnan(values[i])) return not_a_number;
    if (values[i] > largest) largest = values[i];
  }
  return largest;
}

// R's max(x, na.rm = TRUE): an empty or all-missing vector gives -Inf, as R's
// max() does (with a warning R raises and we do not).
double r_max_skip_missing(const double* values, std::int64_t size) {
  double largest = -infinity;
  for (std::int64_t i = 0; i < size; ++i) {
    if (std::isnan(values[i])) continue;
    if (values[i] > largest) largest = values[i];
  }
  return largest;
}

}  // namespace

// ---------------------------------------------------------------------------
// The three full-panel passes.
// ---------------------------------------------------------------------------

Status fit_pass1(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                 Span<const double> solve_x_fixed, std::int64_t p, const CscView& z,
                 const std::vector<SolveBlock>& solve_blocks, Span<const double> beta_in,
                 Span<const double> u_in, Span<const double> lam_diag, const GeneBlock& counts,
                 const GeneBlock& ambient, Span<const double> offset, Span<const double> rho,
                 Span<const double> alpha, Span<const double> sample_weight, bool nb2,
                 bool gaussian, bool seed_iteration, std::int64_t n, std::int64_t q_total,
                 std::int64_t n_genes, std::int64_t chunk_size, std::int64_t sub_genes,
                 int interior_precision, bool last_iter, Span<double> beta_out, Span<double> u_out,
                 Span<double> re_var_out, Span<double> se_beta, Span<double> se_u,
                 std::vector<int>* nan_genes, int n_threads, const InterruptCheck& interrupted) {
  const std::int64_t width = std::min<std::int64_t>(chunk_size, n_genes);
  if (width <= 0) return Status::failure(StatusCode::invalid_argument, "no genes to fit");

  std::vector<double> z_buffer(static_cast<std::size_t>(n * width));
  std::vector<double> w_buffer(static_cast<std::size_t>(n * width));
  std::vector<double> colsum_w(static_cast<std::size_t>(width));
  std::vector<double> eta_scratch;
  std::vector<int> chunk_genes(static_cast<std::size_t>(width));
  std::vector<double> chunk_beta(static_cast<std::size_t>(p * width));
  std::vector<double> chunk_u(static_cast<std::size_t>(q_total * width));
  std::vector<double> chunk_ainv(static_cast<std::size_t>((p + q_total) * width));
  std::vector<double> lam_chunk(static_cast<std::size_t>(q_total * width));

  for (std::int64_t first = 0; first < n_genes; first += width) {
    const std::int64_t m_chunk = std::min(width, n_genes - first);
    for (std::int64_t j = 0; j < m_chunk; ++j) {
      chunk_genes[static_cast<std::size_t>(j)] = static_cast<int>(first + j);
    }

    // ---- the working response, a sub-block of genes at a time ----
    const std::int64_t step = sub_genes > 0 ? std::min<std::int64_t>(sub_genes, m_chunk) : m_chunk;
    for (std::int64_t start = 0; start < m_chunk; start += step) {
      const std::int64_t len = std::min(step, m_chunk - start);
      Span<const double> eta_span;
      if (!seed_iteration) {
        eta_scratch.resize(static_cast<std::size_t>(n * len));
        const Status eta_status =
            eta_block(x1, x1_is_unit, x_fixed, p, beta_in, z, u_in,
                      Span<const int>(chunk_genes.data() + start, len), n, n_genes,
                      Span<double>(eta_scratch.data(), n * len), n_threads, interrupted);
        if (!eta_status.is_ok()) return eta_status;
        eta_span = Span<const double>(eta_scratch.data(), n * len);
      }
      const Status status = working_response(
          eta_span, shifted_block(counts, first + start), shifted_block(ambient, first + start),
          offset, rho, Span<const double>(alpha.data + first + start, len), sample_weight, nb2,
          gaussian, seed_iteration, n, len,
          Span<double>(z_buffer.data() + start * n, n * len),
          Span<double>(w_buffer.data() + start * n, n * len),
          Span<double>(colsum_w.data() + start, len), n_threads, interrupted);
      if (!status.is_ok()) return status;
    }

    // ---- the ridge must not be quantised away in single precision ----
    // Below eps_float * sum(w) the ridge vanishes, at a threshold that falls as
    // 1/n. Decided over the whole logical chunk, as the R it replaces did.
    int chunk_precision = last_iter ? 0 : interior_precision;
    if (chunk_precision != 0) {
      double lam_min = infinity;
      double w_max = -infinity;
      for (std::int64_t j = 0; j < m_chunk; ++j) {
        for (std::int64_t k = 0; k < q_total; ++k) {
          const double value = lam_diag[k + (first + j) * q_total];
          if (value < lam_min) lam_min = value;
        }
        if (colsum_w[static_cast<std::size_t>(j)] > w_max) {
          w_max = colsum_w[static_cast<std::size_t>(j)];
        }
      }
      if (std::isfinite(lam_min) && std::isfinite(w_max) && lam_min < 1e-5 * w_max) {
        chunk_precision = 0;
      }
    }

    // ---- the per-gene solve ----
    for (std::int64_t j = 0; j < m_chunk; ++j) {
      for (std::int64_t k = 0; k < q_total; ++k) {
        lam_chunk[static_cast<std::size_t>(k + j * q_total)] = lam_diag[k + (first + j) * q_total];
      }
    }
    const Status solve_status = solve_genes_chunk(
        solve_x_fixed, n, p, solve_blocks, Span<const double>(w_buffer.data(), n * m_chunk),
        Span<const double>(z_buffer.data(), n * m_chunk),
        Span<const double>(lam_chunk.data(), q_total * m_chunk), q_total, m_chunk,
        chunk_precision != 0, n_threads, interrupted,
        Span<double>(chunk_beta.data(), p * m_chunk),
        Span<double>(chunk_u.data(), q_total * m_chunk),
        Span<double>(chunk_ainv.data(), (p + q_total) * m_chunk));
    if (!solve_status.is_ok()) return solve_status;

    // ---- the lossless NaN guard: re-solve any bad gene in double ----
    // The float per-gene Cholesky returns a NaN BLUP column when a gene's
    // working-weight system is borderline in single precision. One such column
    // makes Z %*% U[, g] non-finite, which NaN-poisons den = rowSums(WA * a) for
    // every cell with a non-zero ambient value and would silently collapse the
    // contamination loading to its prior. The double re-solve is robust because
    // the 1/tau ridge makes every per-gene system positive definite.
    std::vector<std::int64_t> bad;
    for (std::int64_t j = 0; j < m_chunk; ++j) {
      bool finite = true;
      for (std::int64_t k = 0; k < p && finite; ++k) {
        if (!std::isfinite(chunk_beta[static_cast<std::size_t>(k + j * p)])) finite = false;
      }
      for (std::int64_t k = 0; k < q_total && finite; ++k) {
        if (!std::isfinite(chunk_u[static_cast<std::size_t>(k + j * q_total)])) finite = false;
      }
      if (!finite) bad.push_back(j);
    }
    if (!bad.empty() && chunk_precision != 0) {
      const std::int64_t n_bad = static_cast<std::int64_t>(bad.size());
      std::vector<double> w_bad(static_cast<std::size_t>(n * n_bad));
      std::vector<double> z_bad(static_cast<std::size_t>(n * n_bad));
      std::vector<double> lam_bad(static_cast<std::size_t>(q_total * n_bad));
      std::vector<double> beta_bad(static_cast<std::size_t>(p * n_bad));
      std::vector<double> u_bad(static_cast<std::size_t>(q_total * n_bad));
      std::vector<double> ainv_bad(static_cast<std::size_t>((p + q_total) * n_bad));
      for (std::int64_t j = 0; j < n_bad; ++j) {
        const std::int64_t source = bad[static_cast<std::size_t>(j)];
        std::copy(w_buffer.begin() + source * n, w_buffer.begin() + (source + 1) * n,
                  w_bad.begin() + j * n);
        std::copy(z_buffer.begin() + source * n, z_buffer.begin() + (source + 1) * n,
                  z_bad.begin() + j * n);
        std::copy(lam_chunk.begin() + source * q_total, lam_chunk.begin() + (source + 1) * q_total,
                  lam_bad.begin() + j * q_total);
      }
      const Status redo = solve_genes_chunk(
          solve_x_fixed, n, p, solve_blocks, Span<const double>(w_bad.data(), n * n_bad),
          Span<const double>(z_bad.data(), n * n_bad),
          Span<const double>(lam_bad.data(), q_total * n_bad), q_total, n_bad, false, n_threads,
          interrupted, Span<double>(beta_bad.data(), p * n_bad),
          Span<double>(u_bad.data(), q_total * n_bad),
          Span<double>(ainv_bad.data(), (p + q_total) * n_bad));
      if (!redo.is_ok()) return redo;
      for (std::int64_t j = 0; j < n_bad; ++j) {
        const std::int64_t target = bad[static_cast<std::size_t>(j)];
        std::copy(beta_bad.begin() + j * p, beta_bad.begin() + (j + 1) * p,
                  chunk_beta.begin() + target * p);
        std::copy(u_bad.begin() + j * q_total, u_bad.begin() + (j + 1) * q_total,
                  chunk_u.begin() + target * q_total);
        std::copy(ainv_bad.begin() + j * (p + q_total), ainv_bad.begin() + (j + 1) * (p + q_total),
                  chunk_ainv.begin() + target * (p + q_total));
      }
    }

    // ---- write the chunk's columns out ----
    for (std::int64_t j = 0; j < m_chunk; ++j) {
      const std::int64_t gene = first + j;
      bool finite = true;
      for (std::int64_t k = 0; k < p; ++k) {
        const double value = chunk_beta[static_cast<std::size_t>(k + j * p)];
        beta_out[k + gene * p] = value;
        if (!std::isfinite(value)) finite = false;
        if (last_iter) {
          const double variance = chunk_ainv[static_cast<std::size_t>(k + j * (p + q_total))];
          se_beta[k + gene * p] = std::sqrt(variance > 0 ? variance : 0.0);
        }
      }
      for (std::int64_t k = 0; k < q_total; ++k) {
        const double value = chunk_u[static_cast<std::size_t>(k + j * q_total)];
        u_out[k + gene * q_total] = value;
        if (!std::isfinite(value)) finite = false;
        const double variance = chunk_ainv[static_cast<std::size_t>(p + k + j * (p + q_total))];
        re_var_out[k + gene * q_total] = variance > 0 ? variance : 0.0;
        if (last_iter) se_u[k + gene * q_total] = std::sqrt(variance > 0 ? variance : 0.0);
      }
      if (!finite && nan_genes != nullptr) nan_genes->push_back(static_cast<int>(gene + 1));
    }
  }
  return Status::success();
}

Status rho_pass(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                std::int64_t p, const CscView& z, Span<const double> beta_in,
                Span<const double> u_in, Span<const double> prev_beta, Span<const double> prev_u,
                bool have_previous, const GeneBlock& counts, const GeneBlock& ambient,
                Span<const double> offset, Span<const double> rho, Span<const double> alpha,
                Span<const double> mask, Span<const int> mask_index, std::int64_t n_mask_rows,
                bool nb2, std::int64_t n, std::int64_t n_genes, std::int64_t chunk_size,
                Span<double> num, Span<double> den, double* rel_delta_max, double* rel_delta_sum,
                std::int64_t* n_finite, std::int64_t* n_nonfinite, Span<double> tail_counts,
                int n_threads, const InterruptCheck& interrupted) {
  for (std::int64_t i = 0; i < n; ++i) {
    num[i] = 0;
    den[i] = 0;
  }
  *rel_delta_max = 0;
  *rel_delta_sum = 0;
  *n_finite = 0;
  *n_nonfinite = 0;
  for (std::int64_t k = 0; k < tail_counts.size; ++k) tail_counts[k] = 0;

  const std::int64_t width = std::min<std::int64_t>(chunk_size, n_genes);
  if (width <= 0) return Status::failure(StatusCode::invalid_argument, "no genes to accumulate");
  std::vector<double> eta(static_cast<std::size_t>(n * width));
  std::vector<double> previous_eta(have_previous ? static_cast<std::size_t>(n * width) : 0);
  std::vector<int> chunk_genes(static_cast<std::size_t>(width));

  for (std::int64_t first = 0; first < n_genes; first += width) {
    const std::int64_t m_chunk = std::min(width, n_genes - first);
    for (std::int64_t j = 0; j < m_chunk; ++j) {
      chunk_genes[static_cast<std::size_t>(j)] = static_cast<int>(first + j);
    }
    const Span<const int> genes(chunk_genes.data(), m_chunk);
    Status status = eta_block(x1, x1_is_unit, x_fixed, p, beta_in, z, u_in, genes, n, n_genes,
                              Span<double>(eta.data(), n * m_chunk), n_threads, interrupted);
    if (!status.is_ok()) return status;
    if (have_previous) {
      status = eta_block(x1, x1_is_unit, x_fixed, p, prev_beta, z, prev_u, genes, n, n_genes,
                         Span<double>(previous_eta.data(), n * m_chunk), n_threads, interrupted);
      if (!status.is_ok()) return status;
    }
    status = rho_accumulate(
        Span<const double>(eta.data(), n * m_chunk),
        have_previous ? Span<const double>(previous_eta.data(), n * m_chunk) : empty_doubles(),
        shifted_block(counts, first), shifted_block(ambient, first), offset, rho,
        Span<const double>(alpha.data + first, m_chunk), mask, mask_index, n_mask_rows, nb2, false,
        !have_previous, n, m_chunk, num, den, rel_delta_max, rel_delta_sum, n_finite, n_nonfinite,
        tail_counts, n_threads, interrupted);
    if (!status.is_ok()) return status;
  }
  return Status::success();
}

Status dispersion_pass(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                       std::int64_t p, const CscView& z, Span<const double> beta_in,
                       Span<const double> u_in, const GeneBlock& counts, const GeneBlock& ambient,
                       Span<const double> offset, Span<const double> rho, bool nb2, bool gaussian,
                       bool zero_collapse, double max_cells, LogDensity density, bool fast_density,
                       std::int64_t n, std::int64_t n_genes, std::int64_t chunk_size,
                       Span<double> alpha, std::int64_t* n_noninteger, int n_threads,
                       const InterruptCheck& interrupted) {
  *n_noninteger = 0;
  const std::int64_t width = std::min<std::int64_t>(chunk_size, n_genes);
  if (width <= 0) return Status::failure(StatusCode::invalid_argument, "no genes to disperse");
  std::vector<double> eta(static_cast<std::size_t>(n * width));
  std::vector<int> chunk_genes(static_cast<std::size_t>(width));

  for (std::int64_t first = 0; first < n_genes; first += width) {
    const std::int64_t m_chunk = std::min(width, n_genes - first);
    for (std::int64_t j = 0; j < m_chunk; ++j) {
      chunk_genes[static_cast<std::size_t>(j)] = static_cast<int>(first + j);
    }
    const Status eta_status =
        eta_block(x1, x1_is_unit, x_fixed, p, beta_in, z, u_in,
                  Span<const int>(chunk_genes.data(), m_chunk), n, n_genes,
                  Span<double>(eta.data(), n * m_chunk), n_threads, interrupted);
    if (!eta_status.is_ok()) return eta_status;

    if (gaussian) {
      // The identity-link dispersion step: the residual variance, floored where
      // the R engine's pmax(sigma2_new, 1e-8) floors it. No density, no search.
      const Status sigma2_status = residual_variance_chunk(
          Span<const double>(eta.data(), n * m_chunk), shifted_block(counts, first), offset, 1e-8,
          n, m_chunk, Span<double>(alpha.data + first, m_chunk), n_threads, interrupted);
      if (!sigma2_status.is_ok()) return sigma2_status;
      continue;
    }

    std::int64_t chunk_noninteger = 0;
    const Status status = dispersion_chunk(
        Span<const double>(eta.data(), n * m_chunk), shifted_block(counts, first),
        shifted_block(ambient, first), offset, rho, nb2, zero_collapse, max_cells, density,
        fast_density, n, m_chunk, Span<double>(alpha.data + first, m_chunk), &chunk_noninteger,
        n_threads, interrupted);
    if (!status.is_ok()) return status;
    *n_noninteger += chunk_noninteger;
  }
  return Status::success();
}

// ---------------------------------------------------------------------------
// The loop.
// ---------------------------------------------------------------------------

namespace {

// The per-cell contamination fraction the fitting trace prints for the first
// three iterations: the same chunked final pass the fit's own summary makes,
// asked only for the two row sums.
Status contamination_row_sums(const IrlsLoopInputs& inputs, const IrlsLoopOptions& options,
                              const std::vector<double>& beta, const std::vector<double>& u,
                              const std::vector<double>& rho, std::vector<double>& spill_row_sum,
                              std::vector<double>& total_row_sum,
                              const InterruptCheck& interrupted) {
  const std::int64_t n = inputs.n;
  const std::int64_t n_genes = inputs.n_genes;
  const std::int64_t width = std::min<std::int64_t>(options.chunk_size, n_genes);
  spill_row_sum.assign(static_cast<std::size_t>(n), 0.0);
  total_row_sum.assign(static_cast<std::size_t>(n), 0.0);

  // Every cell is in the same nominal group: the diagnostic wants only the row
  // sums, so the per-group block is never asked for.
  const std::vector<int> one_group(static_cast<std::size_t>(n), 0);
  std::vector<double> eta(static_cast<std::size_t>(n * width));
  std::vector<int> chunk_genes(static_cast<std::size_t>(width));
  std::vector<double> mu_column_sum(static_cast<std::size_t>(width));
  std::vector<double> chunk_spill(static_cast<std::size_t>(n));
  std::vector<double> chunk_total(static_cast<std::size_t>(n));

  for (std::int64_t first = 0; first < n_genes; first += width) {
    const std::int64_t m_chunk = std::min(width, n_genes - first);
    for (std::int64_t j = 0; j < m_chunk; ++j) {
      chunk_genes[static_cast<std::size_t>(j)] = static_cast<int>(first + j);
    }
    const Status eta_status = eta_block(
        inputs.x1, inputs.x1_is_unit, inputs.x_fixed, inputs.p,
        Span<const double>(beta.data(), static_cast<std::int64_t>(beta.size())), inputs.z,
        Span<const double>(u.data(), static_cast<std::int64_t>(u.size())),
        Span<const int>(chunk_genes.data(), m_chunk), n, n_genes,
        Span<double>(eta.data(), n * m_chunk), options.n_threads, interrupted);
    if (!eta_status.is_ok()) return eta_status;
    bool any_nonzero = false;
    const Status status = final_pass_statistics(
        Span<const double>(eta.data(), n * m_chunk), n, m_chunk, inputs.offset,
        shifted_block(inputs.ambient, first),
        Span<const double>(rho.data(), static_cast<std::int64_t>(rho.size())),
        Span<const int>(one_group.data(), n), 1, empty_output(), empty_output(),
        Span<double>(mu_column_sum.data(), m_chunk), Span<double>(chunk_spill.data(), n),
        Span<double>(chunk_total.data(), n), empty_output(), empty_output(), &any_nonzero,
        interrupted);
    if (!status.is_ok()) return status;
    for (std::int64_t i = 0; i < n; ++i) {
      spill_row_sum[static_cast<std::size_t>(i)] += chunk_spill[static_cast<std::size_t>(i)];
      total_row_sum[static_cast<std::size_t>(i)] += chunk_total[static_cast<std::size_t>(i)];
    }
  }
  return Status::success();
}

// The variance-component update of one iteration: the EM step, the chosen
// shrinkage, the data-informed prior weights and the cap.
//
// `tau_flat` receives the per-block components the fit reports -- the EM step,
// hierarchically shrunk where that mode asks for it -- and `tau_g_array` the
// per-(column, gene) ridge the next solve inverts. Under the shared and
// hierarchical modes the second is the first broadcast across the genes; under
// the two adaptive modes it is fitted per gene and the two part company, which
// is why both are carried.
Status update_variance_components(const IrlsLoopInputs& inputs, const IrlsLoopOptions& options,
                                  Trigamma trigamma, const std::vector<double>& u,
                                  const std::vector<double>& re_var, std::vector<double>& tau_flat,
                                  std::vector<double>& tau_g_array, std::int64_t* n_binding) {
  const std::int64_t q = inputs.q;
  const std::int64_t n_genes = inputs.n_genes;
  const std::vector<TauBlock>& blocks = *inputs.tau_blocks;

  Status status = tau_em_update(Span<const double>(u.data(), q * n_genes),
                                Span<const double>(re_var.data(), q * n_genes), q, n_genes,
                                Span<double>(tau_flat.data(), q));
  if (!status.is_ok()) return status;

  if (options.tau_shrinkage == TauShrinkage::hierarchical) {
    // R holds the block as a K_terms x K_groups matrix filled BY ROW from the
    // flat slice, so the slice cannot be handed to the shrinkage as it lies:
    // it has to be transposed in and back out.
    for (std::size_t b = 0; b < blocks.size(); ++b) {
      const TauBlock& block = blocks[b];
      const std::int64_t n_terms = block.n_terms;
      const std::int64_t n_groups = block.n_groups;
      std::vector<double> matrix(static_cast<std::size_t>(n_terms * n_groups));
      std::vector<double> shrunk(static_cast<std::size_t>(n_terms * n_groups));
      for (std::int64_t term = 0; term < n_terms; ++term) {
        for (std::int64_t group = 0; group < n_groups; ++group) {
          matrix[static_cast<std::size_t>(term + group * n_terms)] =
              tau_flat[static_cast<std::size_t>(block.col_offset + term * n_groups + group)];
        }
      }
      status = tau_hierarchical(Span<const double>(matrix.data(), n_terms * n_groups),
                                block.group_size, n_terms, n_groups, 0.5,
                                Span<double>(shrunk.data(), n_terms * n_groups));
      if (!status.is_ok()) return status;
      for (std::int64_t term = 0; term < n_terms; ++term) {
        for (std::int64_t group = 0; group < n_groups; ++group) {
          tau_flat[static_cast<std::size_t>(block.col_offset + term * n_groups + group)] =
              shrunk[static_cast<std::size_t>(term + group * n_terms)];
        }
      }
    }
  }

  if (options.tau_shrinkage == TauShrinkage::shared ||
      options.tau_shrinkage == TauShrinkage::hierarchical) {
    for (std::int64_t gene = 0; gene < n_genes; ++gene) {
      for (std::int64_t k = 0; k < q; ++k) {
        tau_g_array[static_cast<std::size_t>(k + gene * q)] = tau_flat[static_cast<std::size_t>(k)];
      }
    }
  } else {
    // s2 = U^2 + V, the per-(column, gene) posterior second moment the two
    // adaptive modes shrink.
    for (std::size_t b = 0; b < blocks.size(); ++b) {
      const TauBlock& block = blocks[b];
      const std::int64_t n_rows = block.n_cols;
      std::vector<double> s2(static_cast<std::size_t>(n_rows * n_genes));
      std::vector<double> shrunk(static_cast<std::size_t>(n_rows * n_genes));
      for (std::int64_t gene = 0; gene < n_genes; ++gene) {
        for (std::int64_t k = 0; k < n_rows; ++k) {
          const std::size_t source = static_cast<std::size_t>(block.col_offset + k + gene * q);
          s2[static_cast<std::size_t>(k + gene * n_rows)] =
              u[source] * u[source] + re_var[source];
        }
      }
      // Wolfinger-O'Connell (1993) / Schall (1991) REML correction on the EM
      // update, off unless R_REML_TAU asks for it.
      double reml_factor = 1;
      if (options.use_reml && block.n_groups > inputs.p) {
        const double denominator = std::max<double>(block.n_groups - inputs.p, 1);
        reml_factor = std::min(block.n_groups / denominator, 2.0);
      }
      if (options.tau_shrinkage == TauShrinkage::adaptive) {
        std::vector<double> panel(static_cast<std::size_t>(n_rows));
        std::vector<double> panel_median(static_cast<std::size_t>(n_rows));
        std::vector<double> log_variance(static_cast<std::size_t>(n_rows));
        std::vector<int> n_log_finite(static_cast<std::size_t>(n_rows));
        status = tau_eb_summaries(Span<const double>(s2.data(), n_rows * n_genes), n_rows, n_genes,
                                  reml_factor, Span<double>(panel.data(), n_rows),
                                  Span<double>(panel_median.data(), n_rows),
                                  Span<double>(log_variance.data(), n_rows),
                                  Span<int>(n_log_finite.data(), n_rows));
        if (!status.is_ok()) return status;
        // Cap the minimum d0: very weak shrinkage produces wild per-gene tau
        // that destabilise the per-gene WLS solve.
        std::vector<double> d0(static_cast<std::size_t>(n_rows));
        for (std::int64_t k = 0; k < n_rows; ++k) {
          const double estimate =
              estimate_d0(log_variance[static_cast<std::size_t>(k)],
                          n_log_finite[static_cast<std::size_t>(k)], 5.0, trigamma);
          d0[static_cast<std::size_t>(k)] = std::max(estimate, options.d0_min);
        }
        status = tau_eb_apply(Span<const double>(s2.data(), n_rows * n_genes), n_rows, n_genes,
                              reml_factor, Span<const double>(panel.data(), n_rows),
                              Span<const double>(panel_median.data(), n_rows),
                              Span<const double>(d0.data(), n_rows), not_a_number,
                              Span<double>(shrunk.data(), n_rows * n_genes));
        if (!status.is_ok()) return status;
      } else {
        status = tau_half_cauchy(Span<const double>(s2.data(), n_rows * n_genes), n_rows, n_genes,
                                 5, 1e-4, empty_doubles(), empty_doubles(),
                                 Span<double>(shrunk.data(), n_rows * n_genes), empty_output(),
                                 empty_output(), empty_output());
        if (!status.is_ok()) return status;
      }
      for (std::int64_t gene = 0; gene < n_genes; ++gene) {
        for (std::int64_t k = 0; k < n_rows; ++k) {
          tau_g_array[static_cast<std::size_t>(block.col_offset + k + gene * q)] =
              shrunk[static_cast<std::size_t>(k + gene * n_rows)];
        }
      }
    }
  }

  if (inputs.data_informed_weights.size != 0) {
    for (std::int64_t i = 0; i < q * n_genes; ++i) {
      tau_g_array[static_cast<std::size_t>(i)] *= inputs.data_informed_weights[i];
    }
    for (std::int64_t i = 0; i < q * n_genes; ++i) {
      tau_g_array[static_cast<std::size_t>(i)] =
          r_pmax(tau_g_array[static_cast<std::size_t>(i)], 1e-8);
    }
  }

  // Every tau write is a pmax against a floor; without a ceiling a column
  // identified only by the ridge climbs by mean(u^2) per iteration, and the
  // inflated posterior variance feeds the next update.
  *n_binding = 0;
  if (std::isfinite(options.tau_max)) {
    status = tau_clamp(Span<double>(tau_g_array.data(), q * n_genes), q * n_genes, options.tau_max,
                       n_binding);
    if (!status.is_ok()) return status;
  }
  return Status::success();
}

}  // namespace

Status run_irls_loop(const IrlsLoopInputs& inputs, const IrlsLoopOptions& options,
                     LogDensity density, Trigamma trigamma, const IrlsLoopReporter& reporter,
                     const IrlsLoopOutputs& outputs, const InterruptCheck& interrupted) {
  if (inputs.solve_blocks == nullptr || inputs.tau_blocks == nullptr) {
    return Status::failure(StatusCode::invalid_argument, "the design blocks are required");
  }
  if (options.n_iter < 1) {
    return Status::failure(StatusCode::invalid_argument, "n_iter must be at least one");
  }
  const std::int64_t n = inputs.n;
  const std::int64_t p = inputs.p;
  const std::int64_t q = inputs.q;
  const std::int64_t n_genes = inputs.n_genes;

  // ---- the state the loop carries between iterations ----
  // Small and full-size only: the coefficients, their posterior variances, the
  // dispersions, the variance components and the contamination loading. No
  // n x G matrix is held across an iteration.
  std::vector<double> beta(static_cast<std::size_t>(p * n_genes), 0.0);
  std::vector<double> u(static_cast<std::size_t>(q * n_genes), 0.0);
  std::vector<double> beta_scratch(static_cast<std::size_t>(p * n_genes), 0.0);
  std::vector<double> u_scratch(static_cast<std::size_t>(q * n_genes), 0.0);
  std::vector<double> re_var(static_cast<std::size_t>(q * n_genes), 0.0);
  std::vector<double> previous_beta(static_cast<std::size_t>(p * n_genes), 0.0);
  std::vector<double> previous_u(static_cast<std::size_t>(q * n_genes), 0.0);
  std::vector<double> alpha(static_cast<std::size_t>(n_genes), 1.0);
  std::vector<double> previous_alpha(static_cast<std::size_t>(n_genes), 1.0);
  std::vector<double> alpha_scratch(static_cast<std::size_t>(n_genes), 0.0);
  std::vector<double> tau_flat(static_cast<std::size_t>(q), 1.0);
  std::vector<double> tau_g_array(static_cast<std::size_t>(q * n_genes), 1.0);
  std::vector<double> lam_diag(static_cast<std::size_t>(q * n_genes), 1.0);
  std::vector<double> rho(static_cast<std::size_t>(n), 0.0);
  std::vector<double> num(static_cast<std::size_t>(n), 0.0);
  std::vector<double> den(static_cast<std::size_t>(n), 0.0);
  std::vector<double> tail_counts(options.rd_diag ? 4 : 0, 0.0);
  std::vector<int> nan_genes;
  if (inputs.alpha_init.size == n_genes) {
    for (std::int64_t j = 0; j < n_genes; ++j) {
      alpha[static_cast<std::size_t>(j)] = inputs.alpha_init[j];
      previous_alpha[static_cast<std::size_t>(j)] = inputs.alpha_init[j];
    }
  }

  std::vector<double> spill_row_sum;
  std::vector<double> total_row_sum;
  std::vector<double> contamination_fraction;

  bool have_previous = false;   // the previous iteration's coefficients exist
  bool converged = false;
  std::int64_t iterations_run = 0;

  for (std::int64_t it = 1; it <= options.n_iter; ++it) {
    const std::chrono::steady_clock::time_point iteration_started =
        std::chrono::steady_clock::now();
    iterations_run = it;
    if (interrupted && interrupted()) {
      return Status::failure(StatusCode::interrupted, "interrupted");
    }

    // ---- the early stop, decided at the TOP from the PREVIOUS iteration ----
    // The lag is deliberate: it makes the stopping iteration the last one, so it
    // runs the standard errors in double and a final dispersion update and the
    // output schema is that of a full run. The MEAN, not the L-infinity max, is
    // the metric: at large n the max is dominated by a handful of jittering
    // cells and never settles even when the bulk has converged. A fit that has
    // lost genes must not stop early and report success, so the previous
    // iteration's non-finite gene count gates it too.
    const bool stop_now =
        it > options.min_iter && options.early_stop_tol > 0 && it >= 2 &&
        std::isfinite(outputs.history_rel_delta_mean[it - 2]) &&
        outputs.history_rel_delta_mean[it - 2] < options.early_stop_tol &&
        outputs.history_n_nan_genes[it - 2] == 0;
    const bool last_iter = (it == options.n_iter) || stop_now;

    for (std::int64_t i = 0; i < q * n_genes; ++i) {
      lam_diag[static_cast<std::size_t>(i)] = 1 / tau_g_array[static_cast<std::size_t>(i)];
    }

    // ---- the gene-chunk solve ----
    nan_genes.clear();
    Status status = fit_pass1(
        inputs.x1, inputs.x1_is_unit, inputs.x_fixed, inputs.solve_x_fixed, p, inputs.z,
        *inputs.solve_blocks, Span<const double>(beta.data(), p * n_genes),
        Span<const double>(u.data(), q * n_genes),
        Span<const double>(lam_diag.data(), q * n_genes), inputs.counts, inputs.ambient,
        inputs.offset, Span<const double>(rho.data(), n),
        Span<const double>(alpha.data(), n_genes), inputs.sample_weight, options.nb2,
        options.gaussian, it == 1, n, q, n_genes, options.chunk_size, options.sub_genes,
        options.interior_precision, last_iter, Span<double>(beta_scratch.data(), p * n_genes),
        Span<double>(u_scratch.data(), q * n_genes), Span<double>(re_var.data(), q * n_genes),
        outputs.se_beta, outputs.se_u, &nan_genes, options.n_threads, interrupted);
    if (!status.is_ok()) return status;
    beta.swap(beta_scratch);
    u.swap(u_scratch);
    if (!nan_genes.empty() && options.verbose && reporter.message) {
      reporter.message(format_text(
          "    [nan-guard] it=%lld: %lld gene(s) STILL non-finite after double solve\n",
          static_cast<long long>(it), static_cast<long long>(nan_genes.size())));
    }

    // ---- the contamination accumulators and the convergence metric ----
    double rel_delta = 0;
    double rd_n = 0;
    double rd_sum = 0;
    double rd_nonfinite = 0;
    double tail_fraction[4] = {0, 0, 0, 0};
    if (!options.gaussian) {
      double pass_max = 0;
      double pass_sum = 0;
      std::int64_t pass_finite = 0;
      std::int64_t pass_nonfinite = 0;
      status = rho_pass(
          inputs.x1, inputs.x1_is_unit, inputs.x_fixed, p, inputs.z,
          Span<const double>(beta.data(), p * n_genes), Span<const double>(u.data(), q * n_genes),
          Span<const double>(have_previous ? previous_beta.data() : beta.data(), p * n_genes),
          Span<const double>(have_previous ? previous_u.data() : u.data(), q * n_genes),
          have_previous, inputs.counts, inputs.ambient, inputs.offset,
          Span<const double>(rho.data(), n), Span<const double>(alpha.data(), n_genes),
          inputs.mask, inputs.mask_index, inputs.n_mask_rows, options.nb2, n, n_genes,
          options.chunk_size, Span<double>(num.data(), n), Span<double>(den.data(), n), &pass_max,
          &pass_sum, &pass_finite, &pass_nonfinite,
          Span<double>(tail_counts.data(), static_cast<std::int64_t>(tail_counts.size())),
          options.n_threads, interrupted);
      if (!status.is_ok()) return status;
      // R's max() propagates a missing value; std::max() would hide it.
      if (std::isnan(pass_max) || pass_max > rel_delta) rel_delta = pass_max;
      rd_n = static_cast<double>(pass_finite);
      rd_sum = pass_sum;
      rd_nonfinite = static_cast<double>(pass_nonfinite);
      for (std::size_t k = 0; k < tail_counts.size(); ++k) tail_fraction[k] = tail_counts[k];
    } else {
      // The Gaussian path never accumulates rho, and its eta-based metric is
      // unusable because the working response does not depend on the fit, so it
      // measures the move in the random-effect coefficients instead. They carry
      // the cell x neighbour spatial signal; the fixed part is intercept-only.
      if (!have_previous) {
        rel_delta = 1e3;
        rd_n = static_cast<double>(q * n_genes);
        rd_sum = 1e3 * static_cast<double>(q * n_genes);
      } else {
        long double sum = 0.0L;
        std::int64_t count = 0;
        double largest = -infinity;
        for (std::int64_t i = 0; i < q * n_genes; ++i) {
          const double previous = previous_u[static_cast<std::size_t>(i)];
          const double change =
              std::fabs(u[static_cast<std::size_t>(i)] - previous) / r_pmax(std::fabs(previous), 1e-3);
          if (!std::isfinite(change)) continue;
          sum += change;
          count += 1;
          if (change > largest) largest = change;
          if (options.rd_diag) {
            if (change > 0.01) tail_fraction[0] += 1;
            if (change > 0.05) tail_fraction[1] += 1;
            if (change > 0.1) tail_fraction[2] += 1;
            if (change > 1.0) tail_fraction[3] += 1;
          }
        }
        rel_delta = largest;
        rd_n = static_cast<double>(count);
        rd_sum = static_cast<double>(sum);
      }
    }

    // ---- the empirical-Bayes shrink of the per-cell contamination loading ----
    if (!options.gaussian) {
      std::int64_t n_nonfinite_cells = 0;
      double rho_bar = 0;
      double den_quantile = 0;
      status = rho_shrink(Span<const double>(num.data(), n), Span<const double>(den.data(), n), n,
                          Span<double>(rho.data(), n), &n_nonfinite_cells, &rho_bar,
                          &den_quantile);
      if (!status.is_ok()) return status;
      if (options.verbose && n_nonfinite_cells > 0 && reporter.message) {
        reporter.message(format_text(
            "    [percell_bleed] it=%lld guarded %lld/%lld non-finite den/num cells "
            "(rho=prior there)\n",
            static_cast<long long>(it), static_cast<long long>(n_nonfinite_cells),
            static_cast<long long>(n)));
      }
      if (options.verbose && it <= 3 && reporter.message) {
        status = contamination_row_sums(inputs, options, beta, u, rho, spill_row_sum,
                                        total_row_sum, interrupted);
        if (!status.is_ok()) return status;
        contamination_fraction.assign(static_cast<std::size_t>(n), 0.0);
        for (std::int64_t i = 0; i < n; ++i) {
          contamination_fraction[static_cast<std::size_t>(i)] =
              spill_row_sum[static_cast<std::size_t>(i)] /
              r_pmax(total_row_sum[static_cast<std::size_t>(i)], 1e-9);
        }
        reporter.message(format_text(
            "    [percell_bleed] it=%lld  rho_i [%.4f,%.4f] med=%.4f  contam_frac med=%.3f "
            "q90=%.3f\n",
            static_cast<long long>(it), r_min(rho.data(), n), r_max(rho.data(), n),
            r_median_strict(rho.data(), n),
            r_median_strict(contamination_fraction.data(), n),
            r_quantile_type7(contamination_fraction.data(), n, 0.9)));
      }
    }

    // ---- the dispersion, unless the warm-up has frozen it ----
    // The NB alpha MLE is about 37% of the per-iteration cost and converges
    // quickly, so after `alpha_warmup` iterations it is frozen until the last
    // one. The Gaussian residual variance is one pass over the data rather than
    // a per-gene Brent search, so there is nothing to amortise and it always
    // updates.
    const bool update_alpha =
        options.gaussian || static_cast<double>(it) <= options.alpha_warmup || last_iter;
    if (update_alpha) {
      std::int64_t n_noninteger = 0;
      status = dispersion_pass(
          inputs.x1, inputs.x1_is_unit, inputs.x_fixed, p, inputs.z,
          Span<const double>(beta.data(), p * n_genes), Span<const double>(u.data(), q * n_genes),
          inputs.counts, inputs.ambient, inputs.offset, Span<const double>(rho.data(), n),
          options.nb2, options.gaussian, options.zero_collapse, options.alpha_max_cells, density,
          options.fast_density, n, n_genes, options.chunk_size,
          Span<double>(alpha_scratch.data(), n_genes), &n_noninteger, options.n_threads,
          interrupted);
      if (!status.is_ok()) return status;
      if (!options.gaussian && n_noninteger > 0 && reporter.warn) {
        reporter.warn(format_text(
            "iter %lld: %.0f gene(s) have counts that are not whole numbers; their dispersion "
            "is undefined and keeps its previous value.",
            static_cast<long long>(it), static_cast<double>(n_noninteger)));
      }
      for (std::int64_t j = 0; j < n_genes; ++j) {
        const double fitted = alpha_scratch[static_cast<std::size_t>(j)];
        alpha[static_cast<std::size_t>(j)] =
            std::isfinite(fitted) ? fitted : previous_alpha[static_cast<std::size_t>(j)];
      }
      // The [1e-4, 50] clamp is the NB dispersion's range. A residual variance
      // is on the scale of the intensities and has no business being clamped to
      // it; the core has already floored it at 1e-8.
      if (!options.gaussian) {
        for (std::int64_t j = 0; j < n_genes; ++j) {
          double value = alpha[static_cast<std::size_t>(j)];
          value = r_pmax(value, 1e-4);
          if (!std::isnan(value) && value > 50) value = 50;
          alpha[static_cast<std::size_t>(j)] = value;
        }
      }
      previous_alpha = alpha;
    } else if (options.verbose && reporter.message) {
      reporter.message(format_text("    [alpha] it=%lld > alpha_warmup=%s: alpha FROZEN\n",
                                   static_cast<long long>(it), options.alpha_warmup_label.c_str()));
    }

    // ---- the variance components ----
    std::int64_t n_binding = 0;
    status = update_variance_components(inputs, options, trigamma, u, re_var, tau_flat,
                                        tau_g_array, &n_binding);
    if (!status.is_ok()) return status;
    if (n_binding > 0 && reporter.warn) {
      reporter.warn(format_text(
          "iter %lld: %lld variance component(s) reached tau_max = %.4g. A binding cap means "
          "the term is identified only by the ridge; inspect the design rather than raising "
          "the cap.",
          static_cast<long long>(it), static_cast<long long>(n_binding), options.tau_max));
    }
    const double tau_max_seen = r_max_skip_missing(tau_g_array.data(), q * n_genes);

    // The next iteration's convergence metric measures the move from here.
    previous_beta = beta;
    previous_u = u;
    have_previous = true;

    // ---- the history ----
    // The early stop reads this back on the next iteration, so it is the only
    // record of the convergence metric the loop keeps: one source of truth for
    // the decision that picks the last iteration.
    for (std::int64_t k = 0; k < q; ++k) {
      outputs.history_tau_flat[k + (it - 1) * q] = tau_flat[static_cast<std::size_t>(k)];
    }
    for (std::int64_t j = 0; j < n_genes; ++j) {
      outputs.history_alpha[j + (it - 1) * n_genes] = alpha[static_cast<std::size_t>(j)];
    }
    const double rel_delta_mean = rd_n > 0 ? rd_sum / rd_n : rel_delta;
    outputs.history_rel_delta[it - 1] = rel_delta;
    outputs.history_rel_delta_mean[it - 1] = rel_delta_mean;
    outputs.history_n_nan_genes[it - 1] = static_cast<double>(nan_genes.size());
    outputs.history_n_nonfinite[it - 1] = rd_nonfinite;
    outputs.history_tau_max_seen[it - 1] = tau_max_seen;
    if (!nan_genes.empty() && reporter.warn_nan_genes) {
      reporter.warn_nan_genes(
          it, Span<const int>(nan_genes.data(), static_cast<std::int64_t>(nan_genes.size())));
    }

    if (options.rd_diag && rd_n > 0 && reporter.message) {
      reporter.message(format_text(
          "  [rd-diag] it=%lld  mean=%.4g  max=%.3g  frac>0.01=%.2e  >0.05=%.2e  >0.1=%.2e  "
          ">1=%.2e  (N=%.2e)\n",
          static_cast<long long>(it), rd_sum / rd_n, rel_delta, tail_fraction[0] / rd_n,
          tail_fraction[1] / rd_n, tail_fraction[2] / rd_n, tail_fraction[3] / rd_n, rd_n));
    }
    if (options.verbose && reporter.message) {
      // Report the max as well as the median: a runaway confined to the
      // intercept row does not move a median taken over the whole block.
      std::string tau_medians;
      for (std::size_t b = 0; b < inputs.tau_blocks->size(); ++b) {
        const TauBlock& block = (*inputs.tau_blocks)[b];
        if (b > 0) tau_medians += ",";
        tau_medians += format_text(
            "%.3f", r_median_strict(tau_flat.data() + block.col_offset, block.n_cols));
      }
      const double elapsed_seconds =
          std::chrono::duration_cast<std::chrono::duration<double> >(
              std::chrono::steady_clock::now() - iteration_started)
              .count();
      reporter.message(format_text(
          "  [mvpql.streaming] iter %lld  rel_delta[mean]=%.3g (max=%.3g)  alpha[med]=%.2f  "
          "tau=[%s] (max %.3g)  (%.1fs)\n",
          static_cast<long long>(it), rel_delta_mean, rel_delta,
          r_median_strict(alpha.data(), n_genes), tau_medians.c_str(), tau_max_seen,
          elapsed_seconds));
    }

    if (stop_now) {
      converged = true;
      if (options.verbose && reporter.message) {
        reporter.message(format_text(
            "  [mvpql.streaming] EARLY STOP at iter %lld (mean rel_delta[%lld]=%.3g < %.3g; "
            "L-inf max=%.3g)\n",
            static_cast<long long>(it), static_cast<long long>(it - 1),
            outputs.history_rel_delta_mean[it - 2], options.early_stop_tol,
            outputs.history_rel_delta[it - 2]));
      }
      break;
    }
  }

  for (std::int64_t i = 0; i < p * n_genes; ++i) outputs.beta[i] = beta[static_cast<std::size_t>(i)];
  for (std::int64_t i = 0; i < q * n_genes; ++i) outputs.u[i] = u[static_cast<std::size_t>(i)];
  for (std::int64_t j = 0; j < n_genes; ++j) outputs.alpha[j] = alpha[static_cast<std::size_t>(j)];
  for (std::int64_t i = 0; i < q * n_genes; ++i) {
    outputs.tau_g_array[i] = tau_g_array[static_cast<std::size_t>(i)];
  }
  for (std::int64_t k = 0; k < q; ++k) outputs.tau_flat[k] = tau_flat[static_cast<std::size_t>(k)];
  for (std::int64_t i = 0; i < n; ++i) outputs.rho[i] = rho[static_cast<std::size_t>(i)];
  if (outputs.iterations_run != nullptr) *outputs.iterations_run = iterations_run;
  if (outputs.converged != nullptr) *outputs.converged = converged;
  return Status::success();
}

}  // namespace pace
