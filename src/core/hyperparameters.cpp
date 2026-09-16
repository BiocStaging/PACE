// hyperparameters.cpp -- implementation of hyperparameters.hpp.
#include "fp_no_contract.hpp"  // must precede the arithmetic below

#include "hyperparameters.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "r_summaries.hpp"

namespace pace {
namespace {

const double kQuietNaN = std::numeric_limits<double>::quiet_NaN();

// vals_g = max(s2_kg, 1e-9) * reml_factor, the vector every tau mode starts from.
void row_values(Span<const double> s2, std::int64_t q, std::int64_t n_genes, std::int64_t k,
                double reml_factor, std::vector<double>& out) {
  out.resize(static_cast<std::size_t>(n_genes));
  for (std::int64_t g = 0; g < n_genes; ++g) {
    out[static_cast<std::size_t>(g)] = r_pmax(s2[k + g * q], 1e-9) * reml_factor;
  }
}

}  // namespace

Status tau_em_update(Span<const double> u, Span<const double> re_var, std::int64_t q,
                     std::int64_t n_genes, Span<double> tau) {
  if (u.size != q * n_genes || re_var.size != q * n_genes || tau.size != q) {
    return Status::failure(StatusCode::invalid_argument, "u and re_var must be q x n_genes");
  }
  std::vector<double> values(static_cast<std::size_t>(n_genes));
  for (std::int64_t k = 0; k < q; ++k) {
    for (std::int64_t g = 0; g < n_genes; ++g) {
      const double blup = u[k + g * q];
      values[static_cast<std::size_t>(g)] = blup * blup + re_var[k + g * q];
    }
    tau[k] = r_pmax(r_mean_array(values.data(), n_genes, true, nullptr), 1e-6);
  }
  return Status::success();
}

Status tau_hierarchical(Span<const double> tau, Span<const int> n_group, std::int64_t n_terms,
                        std::int64_t n_groups, double lambda_factor, Span<double> out) {
  if (tau.size != n_terms * n_groups || out.size != n_terms * n_groups) {
    return Status::failure(StatusCode::invalid_argument, "tau must be n_terms x n_groups");
  }
  if (n_group.size != n_groups) {
    for (std::int64_t k = 0; k < tau.size; ++k) out[k] = tau[k];
    return Status::success();
  }
  std::vector<double> sizes(static_cast<std::size_t>(n_groups));
  for (std::int64_t c = 0; c < n_groups; ++c) {
    sizes[static_cast<std::size_t>(c)] = n_group[c] <= 0 ? 1.0 : static_cast<double>(n_group[c]);
  }
  std::vector<double> scratch(sizes);
  const double lambda = lambda_factor * r_median(scratch);
  std::vector<double> weight(static_cast<std::size_t>(n_groups));
  for (std::int64_t c = 0; c < n_groups; ++c) {
    weight[static_cast<std::size_t>(c)] =
        sizes[static_cast<std::size_t>(c)] / (sizes[static_cast<std::size_t>(c)] + lambda);
  }
  for (std::int64_t t = 0; t < n_terms; ++t) {
    long double weighted = 0.0L;
    long double total = 0.0L;
    for (std::int64_t c = 0; c < n_groups; ++c) {
      const double product = sizes[static_cast<std::size_t>(c)] * tau[t + c * n_terms];
      weighted += product;
      total += sizes[static_cast<std::size_t>(c)];
    }
    const double global = static_cast<double>(weighted) / static_cast<double>(total);
    for (std::int64_t c = 0; c < n_groups; ++c) {
      const double local = weight[static_cast<std::size_t>(c)];
      out[t + c * n_terms] = r_pmax(local * tau[t + c * n_terms] + (1 - local) * global, 1e-6);
    }
  }
  return Status::success();
}

Status tau_eb_summaries(Span<const double> s2, std::int64_t q, std::int64_t n_genes,
                        double reml_factor, Span<double> panel, Span<double> panel_median,
                        Span<double> log_variance, Span<int> n_log_finite) {
  if (s2.size != q * n_genes || panel.size != q || panel_median.size != q ||
      log_variance.size != q || n_log_finite.size != q) {
    return Status::failure(StatusCode::invalid_argument, "sizes do not match q x n_genes");
  }
  std::vector<double> values;
  std::vector<double> logs;
  std::vector<double> scratch;
  for (std::int64_t k = 0; k < q; ++k) {
    row_values(s2, q, n_genes, k, reml_factor, values);
    panel[k] = r_mean_array(values.data(), n_genes, true, nullptr);
    scratch = values;
    panel_median[k] = r_median(scratch);
    logs.clear();
    for (std::int64_t g = 0; g < n_genes; ++g) {
      const double value = std::log(values[static_cast<std::size_t>(g)]);
      if (std::isfinite(value)) logs.push_back(value);
    }
    n_log_finite[k] = static_cast<int>(logs.size());
    log_variance[k] = logs.size() >= 2
                          ? r_variance_array(logs.data(), static_cast<std::int64_t>(logs.size()), false)
                          : kQuietNaN;
  }
  return Status::success();
}

Status tau_eb_apply(Span<const double> s2, std::int64_t q, std::int64_t n_genes,
                    double reml_factor, Span<const double> panel, Span<const double> panel_median,
                    Span<const double> d0, double tau_floor, Span<double> out) {
  if (s2.size != q * n_genes || out.size != q * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "s2 and out must be q x n_genes");
  }
  if (panel.size != q || panel_median.size != q || d0.size != q) {
    return Status::failure(StatusCode::invalid_argument, "per-row inputs must have q entries");
  }
  std::vector<double> values;
  for (std::int64_t k = 0; k < q; ++k) {
    row_values(s2, q, n_genes, k, reml_factor, values);
    const double floor_value =
        std::isfinite(tau_floor) ? tau_floor : r_pmax(panel_median[k] / 100, 1e-4);
    for (std::int64_t g = 0; g < n_genes; ++g) {
      const double shrunk =
          (d0[k] * panel[k] + values[static_cast<std::size_t>(g)]) / (d0[k] + 1);
      out[k + g * q] = r_pmax(shrunk, floor_value);
    }
  }
  return Status::success();
}

Status tau_half_cauchy(Span<const double> s2, std::int64_t q, std::int64_t n_genes, int n_em_iter,
                       double tau_floor, Span<const double> lambda2_prev,
                       Span<const double> a_prev, Span<double> out, Span<double> lambda2_out,
                       Span<double> a_out, Span<double> panel_out) {
  if (s2.size != q * n_genes || out.size != q * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "s2 and out must be q x n_genes");
  }
  const bool warm_lambda = lambda2_prev.size == q;
  const bool warm_a = a_prev.size == q * n_genes;
  std::vector<double> values;
  std::vector<double> scratch;
  std::vector<double> tau(static_cast<std::size_t>(n_genes));
  std::vector<double> a(static_cast<std::size_t>(n_genes));
  std::vector<double> inverse_a(static_cast<std::size_t>(n_genes));
  for (std::int64_t k = 0; k < q; ++k) {
    row_values(s2, q, n_genes, k, 1.0, values);
    const double panel = r_mean_array(values.data(), n_genes, true, nullptr);
    const double floor_value = r_pmax(panel / 100, tau_floor);
    scratch = values;
    double lambda2 = r_pmax(r_median(scratch), 1e-3);
    if (warm_lambda && std::isfinite(lambda2_prev[k])) lambda2 = lambda2_prev[k];
    bool warm_row = warm_a;
    for (std::int64_t g = 0; warm_row && g < n_genes; ++g) {
      if (!std::isfinite(a_prev[k + g * q])) warm_row = false;
    }
    for (std::int64_t g = 0; g < n_genes; ++g) {
      a[static_cast<std::size_t>(g)] = warm_row ? a_prev[k + g * q] : lambda2;
    }
    for (int step = 0; step < n_em_iter; ++step) {
      for (std::int64_t g = 0; g < n_genes; ++g) {
        tau[static_cast<std::size_t>(g)] =
            r_pmax((1 / a[static_cast<std::size_t>(g)] + values[static_cast<std::size_t>(g)] / 2) / 2,
                   floor_value);
      }
      for (std::int64_t g = 0; g < n_genes; ++g) {
        a[static_cast<std::size_t>(g)] =
            r_pmax((1 / tau[static_cast<std::size_t>(g)] + 1 / lambda2) / 2, 1e-9);
        inverse_a[static_cast<std::size_t>(g)] = 1 / a[static_cast<std::size_t>(g)];
      }
      lambda2 = r_pmax(2 * r_mean_array(inverse_a.data(), n_genes, true, nullptr), 1e-4);
    }
    for (std::int64_t g = 0; g < n_genes; ++g) {
      out[k + g * q] = r_pmax(tau[static_cast<std::size_t>(g)], floor_value);
      if (a_out.size != 0) a_out[k + g * q] = a[static_cast<std::size_t>(g)];
    }
    if (lambda2_out.size != 0) lambda2_out[k] = lambda2;
    if (panel_out.size != 0) panel_out[k] = panel;
  }
  return Status::success();
}

Status tau_clamp(Span<double> tau, std::int64_t size, double tau_max, std::int64_t* n_binding) {
  if (tau.size != size) {
    return Status::failure(StatusCode::invalid_argument, "tau size mismatch");
  }
  std::int64_t binding = 0;
  if (!std::isfinite(tau_max)) {
    *n_binding = 0;
    return Status::success();
  }
  for (std::int64_t k = 0; k < size; ++k) {
    if (std::isnan(tau[k])) continue;
    if (tau[k] > tau_max) binding += 1;
  }
  if (binding > 0) {
    for (std::int64_t k = 0; k < size; ++k) {
      if (std::isnan(tau[k])) continue;
      if (tau[k] > tau_max) tau[k] = tau_max;
    }
  }
  *n_binding = binding;
  return Status::success();
}

Status data_informed_weights(Span<const double> detection_rate, std::int64_t n_focals,
                             std::int64_t n_genes, Span<const int> focal_of_row,
                             Span<const double> scale_of_row, std::int64_t q,
                             Span<double> weights) {
  if (detection_rate.size != n_focals * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "detection_rate must be n_focals x n_genes");
  }
  if (focal_of_row.size != q || scale_of_row.size != q || weights.size != q * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "row inputs must have q entries");
  }
  for (std::int64_t k = 0; k < q; ++k) {
    const int focal = focal_of_row[k];
    for (std::int64_t g = 0; g < n_genes; ++g) {
      weights[k + g * q] = focal <= 0 ? 1.0
                                      : detection_rate[(focal - 1) + g * n_focals] * scale_of_row[k];
    }
    // Each row is scaled so its most informative gene has weight one.
    double maximum = -std::numeric_limits<double>::infinity();
    for (std::int64_t g = 0; g < n_genes; ++g) {
      const double value = weights[k + g * q];
      if (std::isnan(value)) continue;
      if (value > maximum) maximum = value;
    }
    if (std::isfinite(maximum) && maximum > 0) {
      for (std::int64_t g = 0; g < n_genes; ++g) weights[k + g * q] = weights[k + g * q] / maximum;
    }
    for (std::int64_t g = 0; g < n_genes; ++g) {
      weights[k + g * q] = r_pmax(weights[k + g * q], 1e-8);
    }
  }
  return Status::success();
}

}  // namespace pace
