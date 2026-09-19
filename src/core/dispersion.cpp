// dispersion.cpp -- implementation of dispersion.hpp.
#include "fp_no_contract.hpp"  // must precede the arithmetic below

#include "dispersion.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <vector>

#include "thread_pool.hpp"

namespace pace {
namespace {

const double kQuietNaN = std::numeric_limits<double>::quiet_NaN();

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

double floor_at(double value, double floor_value) {
  if (std::isnan(value)) return value;
  return value > floor_value ? value : floor_value;
}

}  // namespace

double brent_minimise(double lower, double upper, double tolerance,
                      const std::function<double(double)>& objective,
                      std::int64_t* n_evaluations) {
  // The squared inverse of the golden ratio.
  const double golden = (3. - std::sqrt(5.)) * .5;
  double machine_eps = DBL_EPSILON;
  machine_eps = std::sqrt(machine_eps);

  double a = lower;
  double b = upper;
  double v = a + golden * (b - a);
  double w = v;
  double x = v;
  double d = 0.;
  double e = 0.;
  std::int64_t evaluations = 0;
  auto evaluate = [&](double at) {
    evaluations += 1;
    return objective(at);
  };
  double fx = evaluate(x);
  double fv = fx;
  double fw = fx;
  const double tolerance_third = tolerance / 3.;

  for (;;) {
    const double middle = (a + b) * .5;
    const double tol1 = machine_eps * std::fabs(x) + tolerance_third;
    const double tol2 = tol1 * 2.;
    if (std::fabs(x - middle) <= tol2 - (b - a) * .5) break;

    double p = 0.;
    double q = 0.;
    double r = 0.;
    if (std::fabs(e) > tol1) {  // fit a parabola
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.;
      if (q > 0.) {
        p = -p;
      } else {
        q = -q;
      }
      r = e;
      e = d;
    }

    double u;
    if (std::fabs(p) >= std::fabs(q * .5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
      // a golden-section step
      e = x < middle ? b - x : a - x;
      d = golden * e;
    } else {
      // a parabolic-interpolation step, kept away from the ends
      d = p / q;
      u = x + d;
      if (u - a < tol2 || b - u < tol2) {
        d = tol1;
        if (x >= middle) d = -d;
      }
    }

    if (std::fabs(d) >= tol1) {
      u = x + d;
    } else if (d > 0.) {
      u = x + tol1;
    } else {
      u = x - tol1;
    }
    const double fu = evaluate(u);

    if (fu <= fx) {
      if (u < x) {
        b = x;
      } else {
        a = x;
      }
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    } else {
      if (u < x) {
        a = u;
      } else {
        b = u;
      }
      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  if (n_evaluations != nullptr) *n_evaluations = evaluations;
  return x;
}

double brent_root(double lower, double upper, double tolerance, int max_iterations,
                  const std::function<double(double)>& objective) {
  double a = lower;
  double b = upper;
  double fa = objective(a);
  double fb = objective(b);
  if (fa == 0.0) return a;
  if (fb == 0.0) return b;
  if ((fa > 0) == (fb > 0)) return kQuietNaN;  // no sign change: R raises an error here

  double c = a;
  double fc = fa;
  int iterations = max_iterations;
  while (iterations-- > 0) {
    const double previous_step = b - a;
    if (std::fabs(fc) < std::fabs(fb)) {
      // swap so that b is the best approximation and c the previous one
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    const double tol_act = 2 * DBL_EPSILON * std::fabs(b) + tolerance / 2;
    double new_step = (c - b) / 2;
    if (std::fabs(new_step) <= tol_act || fb == 0.0) return b;

    if (std::fabs(previous_step) >= tol_act && std::fabs(fa) > std::fabs(fb)) {
      // interpolate: linear when only two points are distinct, inverse quadratic otherwise
      double p;
      double q;
      const double cb = c - b;
      const double t1 = fb / fa;
      if (a == c) {
        p = cb * t1;
        q = 1.0 - t1;
      } else {
        const double t2 = fb / fc;
        const double q0 = fa / fc;
        p = t1 * (cb * q0 * (q0 - t2) - (b - a) * (t2 - 1.0));
        q = (q0 - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if (p > 0.0) {
        q = -q;
      } else {
        p = -p;
      }
      if (p < (0.75 * cb * q - std::fabs(tol_act * q) / 2) && p < std::fabs(previous_step * q / 2)) {
        new_step = p / q;
      }
    }
    if (std::fabs(new_step) < tol_act) {
      new_step = new_step > 0 ? tol_act : -tol_act;
    }
    a = b;
    fa = fb;
    b += new_step;
    fb = objective(b);
    if ((fb > 0) == (fc > 0)) {
      c = a;
      fc = fa;
    }
  }
  return b;
}


// The NB1 negative log-likelihood written out, instead of calling a density per
// cell. `Rf_dnbinom_mu` cannot do any of this because it does not know that
// `size` is tied to `mu`; once size = mu/a the sum collapses a long way.
//
//   log f(x; mu/a, mu) = lgamma(x + mu/a) - lgamma(mu/a) - lgamma(x + 1)
//                        - (mu/a) log1p(a) + x (log a - log1p(a))
//
// so summed over cells, with S_mu = sum mu_i and S_x = sum x_i:
//   - the two per-cell logs become the SCALARS log a and log1p(a)
//   - lgamma(x + 1) does not involve a and cannot move the argmin, so it is
//     dropped; this makes the value returned an offset version of the true
//     negative log-likelihood, which is why nothing else may read it
//   - x == 0 gives lgamma(mu/a) - lgamma(mu/a) = 0, the zero_collapse identity
//     the caller already applies
//
// What is left per non-zero cell is lgamma(x + s) - lgamma(s) with s = mu/a.
// For small integer x that is exactly sum_{k<x} log(s + k), which is both
// cheaper than two lgamma calls and better conditioned, since it never forms
// the difference of two large values. Counts on these panels are mostly 0-3.
double nb1_negative_loglik_fast(const double* y, const double* mean, std::int64_t used,
                                double log_alpha) {
  const double a = std::exp(log_alpha);
  const double log1p_a = std::log1p(a);
  const double log_a = log_alpha;          // log(exp(log_alpha)) without the round trip
  const std::int64_t small_count_limit = 8;
  long double total = 0.0L;
  long double count_sum = 0.0L;
  long double mean_sum = 0.0L;
  for (std::int64_t i = 0; i < used; ++i) {
    const double count = y[static_cast<std::size_t>(i)];
    const double value = mean[static_cast<std::size_t>(i)];
    mean_sum += value;
    if (count == 0) continue;              // its two lgamma terms cancel exactly
    count_sum += count;
    const double size = value / a;
    const std::int64_t whole = static_cast<std::int64_t>(count);
    if (static_cast<double>(whole) == count && whole <= small_count_limit) {
      double ratio = 0.0;
      for (std::int64_t k = 0; k < whole; ++k) ratio += std::log(size + static_cast<double>(k));
      total += ratio;
    } else {
      total += std::lgamma(count + size) - std::lgamma(size);
    }
  }
  // The parts that depend on the cells only through their sums. Zero counts
  // need no special case here: their whole contribution IS -(mu/a) log1p(a),
  // which `mean_sum` already carries, so the caller's zero_collapse identity is
  // subsumed rather than applied on top.
  const double sum_mean = static_cast<double>(mean_sum);
  const double sum_count = static_cast<double>(count_sum);
  const double result = static_cast<double>(total) + sum_count * (log_a - log1p_a)
                        - sum_mean * log1p_a / a;
  return -result;
}

double dispersion_mle(Span<const double> counts, Span<const double> mu, bool nb2,
                      bool zero_collapse, double max_cells, LogDensity density,
                      bool fast_density) {
  std::int64_t used = counts.size;
  if (used != mu.size) return kQuietNaN;
  if (used < 10) return kQuietNaN;
  std::vector<double> y(counts.data, counts.data + used);
  std::vector<double> mean(mu.data, mu.data + used);
  // Even subsampling, as seq.int(1, n, length.out = max_cells) selected.
  if (std::isfinite(max_cells) && static_cast<double>(used) > max_cells) {
    const std::int64_t wanted = static_cast<std::int64_t>(max_cells);
    const double step = static_cast<double>(used - 1) / static_cast<double>(wanted - 1);
    std::vector<double> y_keep(static_cast<std::size_t>(wanted));
    std::vector<double> mean_keep(static_cast<std::size_t>(wanted));
    for (std::int64_t k = 0; k < wanted; ++k) {
      const std::int64_t index = static_cast<std::int64_t>(1 + static_cast<double>(k) * step) - 1;
      y_keep[static_cast<std::size_t>(k)] = y[static_cast<std::size_t>(index)];
      mean_keep[static_cast<std::size_t>(k)] = mean[static_cast<std::size_t>(index)];
    }
    y.swap(y_keep);
    mean.swap(mean_keep);
    used = wanted;
  }

  // The zero-count block of the NB1 likelihood in closed form:
  // log f(0; mu/a, mu) = -(mu/a) log1p(a), which depends on the cell only
  // through mu, so the whole block is one term.
  long double zero_mean_sum = 0.0L;
  std::int64_t n_zero = 0;
  if (zero_collapse && !nb2) {
    for (std::int64_t i = 0; i < used; ++i) {
      if (y[static_cast<std::size_t>(i)] == 0) {
        zero_mean_sum += mean[static_cast<std::size_t>(i)];
        n_zero += 1;
      }
    }
  }
  const double mu_zero_sum = static_cast<double>(zero_mean_sum);
  const bool collapse = zero_collapse && !nb2 && n_zero > 0;

  auto negative_log_likelihood = [&](double log_alpha) {
    if (fast_density && !nb2) {
      return nb1_negative_loglik_fast(y.data(), mean.data(), used, log_alpha);
    }
    const double a = std::exp(log_alpha);
    long double total = 0.0L;
    for (std::int64_t i = 0; i < used; ++i) {
      const double count = y[static_cast<std::size_t>(i)];
      if (collapse && count == 0) continue;
      const double value = mean[static_cast<std::size_t>(i)];
      const double size = nb2 ? 1 / a : value / a;
      total += density(count, size, value);
    }
    double result = static_cast<double>(total);
    if (collapse) result = result - mu_zero_sum * std::log1p(a) / a;
    return -result;
  };
  const double tolerance = std::pow(DBL_EPSILON, 0.25);
  return std::exp(brent_minimise(-6.0, 4.0, tolerance, negative_log_likelihood, nullptr));
}

Status dispersion_chunk(Span<const double> eta, const GeneBlock& counts, const GeneBlock& ambient,
                        Span<const double> offset, Span<const double> rho, bool nb2,
                        bool zero_collapse, double max_cells, LogDensity density,
                        bool fast_density, std::int64_t n,
                        std::int64_t n_genes, Span<double> alpha, std::int64_t* n_noninteger,
                        int n_threads, const InterruptCheck& interrupted) {
  if (eta.size != n * n_genes || offset.size != n || rho.size != n || alpha.size != n_genes) {
    return Status::failure(StatusCode::invalid_argument, "dispersion inputs have the wrong size");
  }
  if (density == nullptr) {
    return Status::failure(StatusCode::invalid_argument, "no density function was supplied");
  }
  // One flag per gene, so the count does not depend on the thread count.
  std::vector<char> noninteger(static_cast<std::size_t>(n_genes), 0);

  auto body = [&](std::int64_t begin, std::int64_t end) {
    std::vector<double> y_column;
    std::vector<double> ambient_column;
    std::vector<double> y;
    std::vector<double> mu;
    for (std::int64_t j = begin; j < end; ++j) {
      expand_column(counts, j, n, y_column);
      expand_column(ambient, j, n, ambient_column);
      const double* eta_column = eta.data + j * n;
      y.clear();
      mu.clear();
      bool whole_numbers = true;
      for (std::int64_t i = 0; i < n; ++i) {
        const double mu_bio = floor_at(std::exp(eta_column[i] + offset[i]), 1e-6);
        const double spill = floor_at(ambient_column[static_cast<std::size_t>(i)] * rho[i], 0.0);
        const double fitted = floor_at(mu_bio + spill, 1e-6);
        if (!std::isfinite(fitted) || !(fitted > 1e-8)) continue;
        const double count = y_column[static_cast<std::size_t>(i)];
        if (count != std::floor(count)) whole_numbers = false;
        y.push_back(count);
        mu.push_back(fitted);
      }
      if (!whole_numbers) {
        noninteger[static_cast<std::size_t>(j)] = 1;
        alpha[j] = kQuietNaN;
        continue;
      }
      alpha[j] = dispersion_mle(Span<const double>(y.data(), static_cast<std::int64_t>(y.size())),
                                Span<const double>(mu.data(), static_cast<std::int64_t>(mu.size())),
                                nb2, zero_collapse, max_cells, density, fast_density);
    }
  };
  const Status status = parallel_for(n_genes, n_threads, 1, body, interrupted);
  if (!status.is_ok()) return status;
  *n_noninteger = 0;
  for (std::int64_t j = 0; j < n_genes; ++j) {
    if (noninteger[static_cast<std::size_t>(j)] != 0) *n_noninteger += 1;
  }
  return Status::success();
}

double estimate_d0(double log_variance, std::int64_t n_used, double d0_max, Trigamma trigamma) {
  if (n_used < 5 || !std::isfinite(log_variance)) return d0_max;
  const double excess = log_variance - trigamma(0.5);
  if (excess <= 0) return d0_max;
  const double tolerance = std::pow(DBL_EPSILON, 0.25);
  const double root = brent_root(1e-4, 1e4, tolerance, 1000,
                                 [&](double x) { return trigamma(x) - excess; });
  if (!std::isfinite(root)) return d0_max;
  return std::min(2 * root, d0_max);
}

Status residual_variance_chunk(Span<const double> eta, const GeneBlock& counts,
                               Span<const double> offset, double floor_value, std::int64_t n,
                               std::int64_t n_genes, Span<double> sigma2, int n_threads,
                               const InterruptCheck& interrupted) {
  if (eta.size != n * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "eta must be n * n_genes");
  }
  if (offset.size != n) {
    return Status::failure(StatusCode::invalid_argument, "offset must have one entry per cell");
  }
  if (sigma2.size != n_genes) {
    return Status::failure(StatusCode::invalid_argument, "sigma2 must have one entry per gene");
  }
  for (std::int64_t i = 0; i < n; ++i) {
    if (!std::isfinite(offset[i])) {
      return Status::failure(StatusCode::invalid_argument, "offset must be finite");
    }
  }

  auto body = [&](std::int64_t begin, std::int64_t end) {
    static thread_local std::vector<double> y_column;
    for (std::int64_t j = begin; j < end; ++j) {
      expand_column(counts, j, n, y_column);
      const double* eta_column = eta.data + j * n;
      long double square_sum = 0.0L;
      for (std::int64_t i = 0; i < n; ++i) {
        // Identity link: the fitted mean is eta + offset, with no exp and no
        // lower clip. R forms mu, subtracts, squares, then takes colMeans.
        const double residual = y_column[static_cast<std::size_t>(i)] - (eta_column[i] + offset[i]);
        square_sum += static_cast<long double>(residual * residual);
      }
      const double mean_square = static_cast<double>(square_sum / static_cast<long double>(n));
      sigma2[j] = mean_square > floor_value ? mean_square : floor_value;
    }
  };
  // One gene per block: gene j reads only column j and writes only sigma2[j],
  // and the long double accumulator never crosses a gene boundary, so the
  // answer does not depend on the worker count.
  return parallel_for(n_genes, n_threads, 1, body, interrupted);
}

Status marginal_variance(const GeneBlock& counts, double floor_value, std::int64_t n,
                         std::int64_t n_genes, Span<double> sigma2, int n_threads,
                         const InterruptCheck& interrupted) {
  if (sigma2.size != n_genes) {
    return Status::failure(StatusCode::invalid_argument, "sigma2 must have one entry per gene");
  }
  if (n <= 0) {
    return Status::failure(StatusCode::invalid_argument, "there must be at least one cell");
  }

  auto body = [&](std::int64_t begin, std::int64_t end) {
    static thread_local std::vector<double> y_column;
    for (std::int64_t j = begin; j < end; ++j) {
      expand_column(counts, j, n, y_column);
      long double sum = 0.0L;
      long double square_sum = 0.0L;
      for (std::int64_t i = 0; i < n; ++i) {
        const double y = y_column[static_cast<std::size_t>(i)];
        sum += static_cast<long double>(y);
        square_sum += static_cast<long double>(y * y);
      }
      const long double cells = static_cast<long double>(n);
      const double mean = static_cast<double>(sum / cells);
      const double mean_square = static_cast<double>(square_sum / cells);
      const double variance = mean_square - mean * mean;
      sigma2[j] = variance > floor_value ? variance : floor_value;
    }
  };
  // One gene per block, as residual_variance_chunk partitions it, so the answer
  // does not depend on the worker count.
  return parallel_for(n_genes, n_threads, 1, body, interrupted);
}

}  // namespace pace
