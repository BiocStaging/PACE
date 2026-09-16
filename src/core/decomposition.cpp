// decomposition.cpp -- implementation of decomposition.hpp.
#include "fp_no_contract.hpp"  // must precede the arithmetic below

#include "decomposition.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "r_summaries.hpp"
#include "thread_pool.hpp"

namespace pace {
namespace {

const double kQuietNaN = std::numeric_limits<double>::quiet_NaN();

// sum_k a_k b_k as R's colSums() of the elementwise product: the product in
// double, the sum in long double. `skip_missing` reproduces na.rm = TRUE.
double dot_product(const double* a, const double* b, std::int64_t n, bool skip_missing) {
  long double sum = 0.0L;
  for (std::int64_t k = 0; k < n; ++k) {
    const double product = a[k] * b[k];
    if (skip_missing && std::isnan(product)) continue;
    sum += product;
  }
  return static_cast<double>(sum);
}

// The quadratic form b' S b for one gene, as R computes
// colSums(B * (S %*% B)) + colSums(SE2 * diag(S), na.rm = TRUE):
// S %*% B in double (a matrix product), the outer sums in long double.
double quadratic_form(const double* s, const double* b, const double* se_squared, std::int64_t p,
                      std::vector<double>& work) {
  work.resize(static_cast<std::size_t>(p));
  for (std::int64_t k = 0; k < p; ++k) {
    double sum = 0.0;
    for (std::int64_t l = 0; l < p; ++l) sum += s[k + l * p] * b[l];
    work[static_cast<std::size_t>(k)] = sum;
  }
  const double centre = dot_product(b, work.data(), p, false);
  if (se_squared == nullptr) return centre;
  long double noise = 0.0L;
  for (std::int64_t k = 0; k < p; ++k) {
    const double product = se_squared[k] * s[k + k * p];
    if (std::isnan(product)) continue;
    noise += product;
  }
  return centre + static_cast<double>(noise);
}

// The values of one 1-based column of a sparse matrix at the 1-based `rows`, which
// must be strictly increasing (they come from a column's own row indices). Both
// sides are walked once, so nothing of the matrix's height is allocated.
Status gather_column(const CscView& matrix, Span<const int> rows, int column, Span<double> values) {
  // checked before the shift, so R's NA_integer_ cannot overflow on the way in
  if (column < 1 || column > matrix.n_cols) {
    return Status::failure(StatusCode::invalid_argument, "column index out of range");
  }
  const int j = column - 1;
  for (std::int64_t r = 0; r < rows.size; ++r) {
    if (rows[r] < 1 || rows[r] > matrix.n_rows) {
      return Status::failure(StatusCode::invalid_argument, "row index out of range");
    }
    if (r > 0 && rows[r] <= rows[r - 1]) {
      return Status::failure(StatusCode::invalid_argument, "rows must be strictly increasing");
    }
    values[r] = 0.0;
  }
  std::int64_t r = 0;
  for (int k = matrix.column_pointer[j]; k < matrix.column_pointer[j + 1] && r < rows.size; ++k) {
    const int row = matrix.row_index[k] + 1;
    while (r < rows.size && rows[r] < row) ++r;
    if (r < rows.size && rows[r] == row) {
      values[r] = matrix.values[k];
      ++r;
    }
  }
  return Status::success();
}

}  // namespace

Status variance_decomposition(const DecompositionInput& in, const DecompositionOutput& out,
                              int n_threads, const InterruptCheck& interrupted) {
  const std::int64_t n_genes = in.n_genes;
  const std::int64_t n_focals = in.n_focals;
  const std::int64_t n_groups = in.n_groups;
  const std::int64_t rows = n_focals * n_genes;
  if (in.ct_means.size != n_groups * n_genes || in.group_size.size != n_groups) {
    return Status::failure(StatusCode::invalid_argument, "ct_means must be n_groups x n_genes");
  }
  if (in.focal_group.size != n_focals || in.n_focal.size != n_focals ||
      in.intercept_rows.size != n_focals) {
    return Status::failure(StatusCode::invalid_argument, "per-focal inputs must have n_focals entries");
  }
  if (in.u.size != in.q * n_genes || in.se_u.size != in.q * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "u and se_u must be q x n_genes");
  }
  if (in.slope_rows.size != in.n_kernel * n_focals ||
      in.kernel_cov.size != n_groups * in.n_kernel * in.n_kernel) {
    return Status::failure(StatusCode::invalid_argument, "kernel slope rows or covariances have the wrong size");
  }
  if (in.n_responder > 0 &&
      (in.responder_rows.size != in.n_responder * n_focals || in.responder_keep.size != in.n_responder ||
       in.responder_cov.size != n_groups * in.n_kernel * in.n_kernel)) {
    return Status::failure(StatusCode::invalid_argument, "responder rows or covariances have the wrong size");
  }
  if (in.n_spill > 0 &&
      (in.spill_cov.size != n_groups * in.n_spill * in.n_spill || in.beta_spill.size != in.n_spill * n_genes)) {
    return Status::failure(StatusCode::invalid_argument, "spillover covariances have the wrong size");
  }
  if (in.toff_var.size != 0 && in.toff_var.size != rows) {
    return Status::failure(StatusCode::invalid_argument, "toff_var must be n_focals x n_genes");
  }
  if (in.mu_mean.size != rows || in.alpha.size != n_genes) {
    return Status::failure(StatusCode::invalid_argument, "mu_mean or alpha have the wrong size");
  }
  if (out.v_state_baseline.size != rows || out.keep.size != rows) {
    return Status::failure(StatusCode::invalid_argument, "outputs must be n_focals x n_genes");
  }
  // Every index the loop below dereferences, checked once here.
  for (std::int64_t f = 0; f < n_focals; ++f) {
    if (in.focal_group[f] < 1 || in.focal_group[f] > n_groups) {
      return Status::failure(StatusCode::invalid_argument, "focal_group index out of range");
    }
    if (in.intercept_rows[f] < 1 || in.intercept_rows[f] > in.q) {
      return Status::failure(StatusCode::invalid_argument, "intercept row index out of range");
    }
    for (std::int64_t k = 0; k < in.n_kernel; ++k) {
      const int row = in.slope_rows[k + f * in.n_kernel];
      if (row < 1 || row > in.q) {
        return Status::failure(StatusCode::invalid_argument, "slope row index out of range");
      }
    }
    for (std::int64_t a = 0; a < in.n_responder; ++a) {
      const int row = in.responder_rows[a + f * in.n_responder];
      if (row < 1 || row > in.q) {
        return Status::failure(StatusCode::invalid_argument, "responder row index out of range");
      }
    }
  }
  for (std::int64_t a = 0; a < in.n_responder; ++a) {
    if (in.responder_keep[a] < 1 || in.responder_keep[a] > in.n_kernel) {
      return Status::failure(StatusCode::invalid_argument, "responder keep index out of range");
    }
  }

  // A cell type with no cells contributes a zero mean count, as the R code set.
  auto type_mean = [&](std::int64_t group, std::int64_t gene) {
    if (in.group_size[group] == 0) return 0.0;
    return in.ct_means[group + gene * n_groups];
  };

  auto body = [&](std::int64_t begin, std::int64_t end) {
    std::vector<double> slope(static_cast<std::size_t>(in.n_kernel));
    std::vector<double> slope_se_squared(static_cast<std::size_t>(in.n_kernel));
    std::vector<double> responder_slope(static_cast<std::size_t>(std::max<std::int64_t>(in.n_responder, 1)));
    std::vector<double> responder_se_squared(responder_slope.size());
    std::vector<double> responder_cov_block(
        static_cast<std::size_t>(std::max<std::int64_t>(in.n_responder * in.n_responder, 1)));
    std::vector<double> spill_slope(static_cast<std::size_t>(std::max<std::int64_t>(in.n_spill, 1)));
    std::vector<double> work;
    for (std::int64_t gene = begin; gene < end; ++gene) {
      for (std::int64_t f = 0; f < n_focals; ++f) {
        const std::int64_t slot = f * n_genes + gene;
        const std::int64_t group = in.focal_group[f] - 1;
        if (in.n_focal[f] < 5 || group < 0) {
          out.v_state_baseline[slot] = kQuietNaN;
          out.v_state_responder[slot] = kQuietNaN;
          out.v_spill[slot] = kQuietNaN;
          out.v_disp[slot] = kQuietNaN;
          out.celltype_offset_sq[slot] = kQuietNaN;
          out.focal_mean[slot] = kQuietNaN;
          out.max_other_mean[slot] = kQuietNaN;
          out.spec[slot] = kQuietNaN;
          out.focal_other_ratio[slot] = kQuietNaN;
          out.is_contaminated[slot] = 0;
          out.pct_celltype[slot] = kQuietNaN;
          out.pct_state[slot] = kQuietNaN;
          out.pct_responder[slot] = kQuietNaN;
          out.pct_spill[slot] = kQuietNaN;
          out.pct_residual[slot] = kQuietNaN;
          out.keep[slot] = 0;
          continue;
        }
        // Specificity against the loudest other cell type (max over types, na.rm).
        const double focal_mean = type_mean(group, gene);
        double other_max = -std::numeric_limits<double>::infinity();
        for (std::int64_t g = 0; g < n_groups; ++g) {
          if (g == group) continue;
          const double value = type_mean(g, gene);
          if (std::isnan(value)) continue;
          if (value > other_max) other_max = value;
        }
        const double spec = focal_mean / r_pmax(focal_mean + other_max, 1e-9);
        const double focal_other_ratio = focal_mean / r_pmax(other_max, 1e-9);

        // Baseline spatial state: a quadratic form in the slope BLUPs.
        for (std::int64_t k = 0; k < in.n_kernel; ++k) {
          const std::int64_t row = in.slope_rows[k + f * in.n_kernel] - 1;
          const double se = in.se_u[row + gene * in.q];
          slope[static_cast<std::size_t>(k)] = in.u[row + gene * in.q];
          slope_se_squared[static_cast<std::size_t>(k)] = se * se;
        }
        const double* kernel_cov_block = in.kernel_cov.data + group * in.n_kernel * in.n_kernel;
        const double v_state_baseline = quadratic_form(kernel_cov_block, slope.data(),
                                                       slope_se_squared.data(), in.n_kernel, work);

        // Responder spatial state: the same form over the condition interaction terms.
        double v_state_responder = 0.0;
        if (in.n_responder > 0) {
          const double* responder_full = in.responder_cov.data + group * in.n_kernel * in.n_kernel;
          for (std::int64_t a = 0; a < in.n_responder; ++a) {
            const std::int64_t row = in.responder_rows[a + f * in.n_responder] - 1;
            const double se = in.se_u[row + gene * in.q];
            responder_slope[static_cast<std::size_t>(a)] = in.u[row + gene * in.q];
            responder_se_squared[static_cast<std::size_t>(a)] = se * se;
            for (std::int64_t b = 0; b < in.n_responder; ++b) {
              responder_cov_block[static_cast<std::size_t>(a + b * in.n_responder)] =
                  responder_full[(in.responder_keep[a] - 1) + (in.responder_keep[b] - 1) * in.n_kernel];
            }
          }
          v_state_responder = quadratic_form(responder_cov_block.data(), responder_slope.data(),
                                             responder_se_squared.data(), in.n_responder, work);
        }

        // Spillover: the stored offset variance, else the legacy _near effects.
        double v_spill = 0.0;
        if (in.toff_var.size != 0) {
          v_spill = in.toff_var[f + gene * n_focals];
        } else if (in.n_spill > 0) {
          for (std::int64_t k = 0; k < in.n_spill; ++k) {
            spill_slope[static_cast<std::size_t>(k)] = in.beta_spill[k + gene * in.n_spill];
          }
          const double* spill_block = in.spill_cov.data + group * in.n_spill * in.n_spill;
          v_spill = quadratic_form(spill_block, spill_slope.data(), nullptr, in.n_spill, work);
        }

        // Dispersion floor of the fitted mean (NB1 is canonical; Leckie 2020).
        const double mu_bar = r_pmax(in.mu_mean[f + gene * n_focals], 1e-9);
        const double alpha = r_pmax(in.alpha[gene], 0.0);
        // R's log(1 + x), not log1p(x): the two differ in the last bits.
        const double v_disp = in.nb1 ? std::log(1 + (1 + alpha) / mu_bar)
                                     : std::log(1 + 1 / mu_bar + alpha);

        const std::int64_t intercept_row = in.intercept_rows[f] - 1;
        const double intercept_u = in.u[intercept_row + gene * in.q];
        const double intercept_se = in.se_u[intercept_row + gene * in.q];
        const double celltype_offset_sq = intercept_u * intercept_u + intercept_se * intercept_se;

        const double total = celltype_offset_sq + v_state_baseline + v_state_responder + v_spill + v_disp;
        out.v_state_baseline[slot] = v_state_baseline;
        out.v_state_responder[slot] = v_state_responder;
        out.v_spill[slot] = v_spill;
        out.v_disp[slot] = v_disp;
        out.celltype_offset_sq[slot] = celltype_offset_sq;
        out.focal_mean[slot] = focal_mean;
        out.max_other_mean[slot] = other_max;
        out.spec[slot] = std::isfinite(spec) ? spec : 0.0;
        out.focal_other_ratio[slot] = focal_other_ratio;
        out.is_contaminated[slot] = (!std::isfinite(focal_other_ratio) || focal_other_ratio < 1) ? 1 : 0;
        out.pct_celltype[slot] = celltype_offset_sq / total * 100;
        out.pct_state[slot] = v_state_baseline / total * 100;
        out.pct_responder[slot] = v_state_responder / total * 100;
        out.pct_spill[slot] = v_spill / total * 100;
        out.pct_residual[slot] = v_disp / total * 100;
        out.keep[slot] = (std::isfinite(total) && total > 0) ? 1 : 0;
      }
    }
  };
  return parallel_for(n_genes, n_threads, 8, body, interrupted);
}

Status decomposition_aggregates(Span<const int> focal_code, int n_focals, std::int64_t n_rows,
                                Span<const double> pct, Span<const double> components,
                                Span<const double> spec, Span<const int> is_contaminated,
                                Span<double> mean_table, Span<double> specw_table,
                                Span<double> pooled_table, Span<double> total,
                                Span<double> total_ss, Span<int> n_rows_out, Span<int> n_specific) {
  const int n_blocks = 5;
  if (focal_code.size != n_rows || pct.size != n_rows * n_blocks || components.size != n_rows * n_blocks ||
      spec.size != n_rows || is_contaminated.size != n_rows || total.size != n_rows) {
    return Status::failure(StatusCode::invalid_argument, "row inputs must all have n_rows entries");
  }
  if (mean_table.size != static_cast<std::int64_t>(n_focals) * n_blocks ||
      pooled_table.size != static_cast<std::int64_t>(n_focals) * n_blocks ||
      (specw_table.size != 0 && specw_table.size != static_cast<std::int64_t>(n_focals) * n_blocks)) {
    return Status::failure(StatusCode::invalid_argument, "tables must be n_focals x 5");
  }
  if (total_ss.size != n_focals || n_rows_out.size != n_focals || n_specific.size != n_focals) {
    return Status::failure(StatusCode::invalid_argument, "per-focal outputs must have n_focals entries");
  }

  // Total = V_state_baseline + V_state_responder + V_spill + V_disp + celltype_offset_sq,
  // in the order the R mutate() added them.
  for (std::int64_t r = 0; r < n_rows; ++r) {
    total[r] = components[r + 1 * n_rows] + components[r + 2 * n_rows] + components[r + 3 * n_rows] +
               components[r + 4 * n_rows] + components[r + 0 * n_rows];
  }

  std::vector<std::vector<std::int64_t>> rows_of_focal(n_focals);
  for (std::int64_t r = 0; r < n_rows; ++r) {
    const int f = focal_code[r];
    if (f >= 0 && f < n_focals) rows_of_focal[f].push_back(r);
  }
  std::vector<double> spec_squared(static_cast<std::size_t>(n_rows));
  std::vector<double> weighted(static_cast<std::size_t>(n_rows));
  for (std::int64_t r = 0; r < n_rows; ++r) spec_squared[static_cast<std::size_t>(r)] = spec[r] * spec[r];

  for (int f = 0; f < n_focals; ++f) {
    const std::vector<std::int64_t>& rows = rows_of_focal[f];
    n_rows_out[f] = static_cast<int>(rows.size());
    int specific = 0;
    for (std::int64_t r : rows) {
      if (is_contaminated[r] == 0) specific += 1;
    }
    n_specific[f] = specific;
    const double spec_weight = static_cast<double>(r_sum(spec_squared.data(), rows, true));
    const double pooled_total = static_cast<double>(r_sum(total.data, rows, true));
    total_ss[f] = pooled_total;
    for (int b = 0; b < n_blocks; ++b) {
      const double* pct_column = pct.data + b * n_rows;
      const double* component_column = components.data + b * n_rows;
      mean_table[f + b * n_focals] = r_mean(pct_column, rows, true, nullptr);
      if (specw_table.size != 0) {
        for (std::int64_t r : rows) {
          weighted[static_cast<std::size_t>(r)] = pct_column[r] * spec_squared[static_cast<std::size_t>(r)];
        }
        specw_table[f + b * n_focals] =
            static_cast<double>(r_sum(weighted.data(), rows, true)) / r_pmax(spec_weight, 1e-12);
      }
      pooled_table[f + b * n_focals] =
          100 * static_cast<double>(r_sum(component_column, rows, true)) / r_pmax(pooled_total, 1e-12);
    }
  }
  return Status::success();
}

Status four_block_shares(Span<const double> celltype_offset_sq, Span<const double> v_state_baseline,
                         Span<const double> v_state_responder, Span<const double> v_spill,
                         Span<const double> v_disp, std::int64_t n_rows, Span<double> total,
                         Span<double> v_state, Span<double> pct_celltype, Span<double> pct_state,
                         Span<double> pct_spill, Span<double> pct_residual) {
  if (celltype_offset_sq.size != n_rows || v_state_baseline.size != n_rows ||
      v_state_responder.size != n_rows || v_spill.size != n_rows || v_disp.size != n_rows ||
      total.size != n_rows || v_state.size != n_rows) {
    return Status::failure(StatusCode::invalid_argument, "every span must have n_rows entries");
  }
  for (std::int64_t r = 0; r < n_rows; ++r) {
    const double sum = celltype_offset_sq[r] + v_state_baseline[r] + v_state_responder[r] + v_spill[r] +
                       v_disp[r];
    const double floor_sum = r_pmax(sum, 1e-12);
    total[r] = sum;
    v_state[r] = v_state_baseline[r] + v_state_responder[r];
    pct_celltype[r] = 100 * celltype_offset_sq[r] / floor_sum;
    pct_state[r] = 100 * v_state[r] / floor_sum;
    pct_spill[r] = 100 * v_spill[r] / floor_sum;
    pct_residual[r] = 100 * v_disp[r] / floor_sum;
  }
  return Status::success();
}

Status single_frame_shares(Span<const double> focal_mean, Span<const double> within_ss,
                           Span<const double> global_mean, Span<const int> group_size, int n_groups,
                           std::int64_t n_genes, Span<const double> v_state,
                           Span<const double> v_responder, Span<const double> v_spill,
                           Span<const double> v_disp, Span<double> pct_celltype,
                           Span<double> pct_spatial, Span<double> pct_responder,
                           Span<double> pct_spill, Span<double> pct_residual, Span<double> ss_lineage,
                           Span<double> ss_within, Span<double> denom, Span<int> keep, int n_threads,
                           const InterruptCheck& interrupted) {
  const std::int64_t slots = static_cast<std::int64_t>(n_groups) * n_genes;
  if (focal_mean.size != slots || within_ss.size != slots || global_mean.size != n_genes ||
      group_size.size != n_groups) {
    return Status::failure(StatusCode::invalid_argument, "statistics must be n_groups x n_genes");
  }
  if (v_state.size != slots || v_spill.size != slots || v_disp.size != slots ||
      (v_responder.size != 0 && v_responder.size != slots)) {
    return Status::failure(StatusCode::invalid_argument, "fit components must be n_groups x n_genes");
  }
  if (keep.size != n_groups) {
    return Status::failure(StatusCode::invalid_argument, "keep must have one entry per group");
  }
  for (int g = 0; g < n_groups; ++g) keep[g] = group_size[g] >= 5 ? 1 : 0;

  auto body = [&](std::int64_t begin, std::int64_t end) {
    for (std::int64_t gene = begin; gene < end; ++gene) {
      for (int g = 0; g < n_groups; ++g) {
        const std::int64_t slot = g + gene * n_groups;
        const std::int64_t row = static_cast<std::int64_t>(g) * n_genes + gene;
        if (keep[g] == 0) {
          pct_celltype[row] = kQuietNaN;
          pct_spatial[row] = kQuietNaN;
          if (pct_responder.size != 0) pct_responder[row] = kQuietNaN;
          pct_spill[row] = kQuietNaN;
          pct_residual[row] = kQuietNaN;
          ss_lineage[row] = kQuietNaN;
          ss_within[row] = kQuietNaN;
          denom[row] = kQuietNaN;
          continue;
        }
        const double deviation = focal_mean[slot] - global_mean[gene];
        const double lineage = group_size[g] * (deviation * deviation);
        const double within = within_ss[slot];
        const double total = lineage + within;
        const double state = v_state[slot];
        const double responder = v_responder.size != 0 ? v_responder[slot] : 0.0;
        const double spill = v_spill[slot];
        const double disp = v_disp[slot];
        const double fit_total = state + responder + spill + disp;
        // ifelse(vt > 0, V / vt, 0): a missing total keeps the row missing.
        const bool missing = std::isnan(fit_total);
        const double p_state = missing ? kQuietNaN : (fit_total > 0 ? state / fit_total : 0.0);
        const double p_responder = missing ? kQuietNaN : (fit_total > 0 ? responder / fit_total : 0.0);
        const double p_spill = missing ? kQuietNaN : (fit_total > 0 ? spill / fit_total : 0.0);
        const double p_disp = missing ? kQuietNaN : (fit_total > 0 ? disp / fit_total : 0.0);
        pct_celltype[row] = 100 * lineage / total;
        pct_spatial[row] = 100 * within * p_state / total;
        if (pct_responder.size != 0) pct_responder[row] = 100 * within * p_responder / total;
        pct_spill[row] = 100 * within * p_spill / total;
        pct_residual[row] = 100 * within * p_disp / total;
        ss_lineage[row] = lineage;
        ss_within[row] = within;
        denom[row] = total;
      }
    }
  };
  return parallel_for(n_genes, n_threads, 16, body, interrupted);
}

Status single_frame_focal_blocks(Span<const int> focal_code, int n_focals, std::int64_t n_rows,
                                 Span<const double> pct, int n_blocks, Span<const double> denom,
                                 Span<double> table) {
  if (focal_code.size != n_rows || pct.size != n_rows * n_blocks || denom.size != n_rows) {
    return Status::failure(StatusCode::invalid_argument, "row inputs must all have n_rows entries");
  }
  if (table.size != static_cast<std::int64_t>(n_focals) * n_blocks) {
    return Status::failure(StatusCode::invalid_argument, "table must be n_focals x n_blocks");
  }
  std::vector<std::vector<std::int64_t>> rows_of_focal(n_focals);
  for (std::int64_t r = 0; r < n_rows; ++r) {
    const int f = focal_code[r];
    if (f >= 0 && f < n_focals) rows_of_focal[f].push_back(r);
  }
  std::vector<double> weighted(static_cast<std::size_t>(n_rows));
  for (int f = 0; f < n_focals; ++f) {
    const std::vector<std::int64_t>& rows = rows_of_focal[f];
    const double total = static_cast<double>(r_sum(denom.data, rows, false));
    for (int b = 0; b < n_blocks; ++b) {
      const double* column = pct.data + b * n_rows;
      for (std::int64_t r : rows) {
        weighted[static_cast<std::size_t>(r)] = column[r] / 100 * denom[r];
      }
      table[f + b * n_focals] = 100 * static_cast<double>(r_sum(weighted.data(), rows, false)) / total;
    }
  }
  return Status::success();
}

Status driver_scores(Span<const double> estimate_shrunk, Span<const double> spec,
                     Span<const double> focal_mean, Span<const double> mu_bar,
                     Span<const double> alpha, Span<const double> u_raw, double var_n, double var_rn,
                     bool has_responder, std::int64_t n, Span<double> mcsd, Span<double> mcsd4,
                     Span<double> v_resid, Span<double> v_s, Span<double> v_rxs, Span<double> v_total,
                     Span<double> r2_s, Span<double> r2_rxs) {
  if (estimate_shrunk.size != n || spec.size != n || focal_mean.size != n || mu_bar.size != n ||
      alpha.size != n) {
    return Status::failure(StatusCode::invalid_argument, "every gene input must have n entries");
  }
  if (has_responder && u_raw.size != n) {
    return Status::failure(StatusCode::invalid_argument, "a condition fit needs the raw BLUPs");
  }
  for (std::int64_t g = 0; g < n; ++g) {
    const double b = estimate_shrunk[g];
    const double b_squared = b * b;
    const double spec_squared = spec[g] * spec[g];
    const double level = r_pmax(focal_mean[g], 0.0);
    mcsd[g] = b_squared * spec_squared * level;
    // spec^4 is R's ^ on a double exponent, which is libm's pow().
    mcsd4[g] = b_squared * std::pow(spec[g], 4.0) * level;
    const double residual = std::log(1 + (1 + r_pmax(alpha[g], 0.0)) / r_pmax(mu_bar[g], 1e-6));
    v_resid[g] = residual;
    if (!has_responder) {
      const double state = b_squared * var_n;
      const double total = state + residual;
      v_s[g] = state;
      v_total[g] = total;
      r2_s[g] = state / r_pmax(total, 1e-12);
      continue;
    }
    const double raw = u_raw[g];
    const double state = raw * raw * var_n;
    const double interaction = b_squared * var_rn;
    const double total = state + interaction + residual;
    v_s[g] = state;
    v_rxs[g] = interaction;
    v_total[g] = total;
    r2_s[g] = state / r_pmax(total, 1e-12);
    r2_rxs[g] = interaction / r_pmax(total, 1e-12);
  }
  return Status::success();
}

Status subset_covariance(const CscView& matrix, Span<const int> rows, Span<const int> cols,
                         Span<const double> scale, Span<double> covariance) {
  const std::int64_t n_obs = rows.size;
  const std::int64_t p = cols.size;
  if (covariance.size != p * p) {
    return Status::failure(StatusCode::invalid_argument, "covariance must be n_cols x n_cols");
  }
  if (scale.size != 0 && scale.size != n_obs) {
    return Status::failure(StatusCode::invalid_argument, "scale must be empty or one value per row");
  }
  if (n_obs < 2) {
    for (std::int64_t k = 0; k < p * p; ++k) covariance[k] = kQuietNaN;
    return Status::success();
  }
  // Dense block of the selected rows and columns, gathered from the sparse matrix
  // by walking each column's stored entries alongside the ascending row list.
  std::vector<double> block(static_cast<std::size_t>(n_obs * p), 0.0);
  for (std::int64_t j = 0; j < p; ++j) {
    const Status gathered = gather_column(matrix, rows, cols[j],
                                          Span<double>(block.data() + j * n_obs, n_obs));
    if (!gathered.is_ok()) return gathered;
  }
  if (scale.size != 0) {
    for (std::int64_t j = 0; j < p; ++j) {
      for (std::int64_t r = 0; r < n_obs; ++r) {
        block[static_cast<std::size_t>(r + j * n_obs)] *= scale[r];
      }
    }
  }
  for (std::int64_t k = 0; k < n_obs * p; ++k) {
    if (!std::isfinite(block[static_cast<std::size_t>(k)])) {
      return Status::failure(StatusCode::invalid_argument, "covariance input must be finite");
    }
  }
  std::vector<std::int64_t> all_rows(static_cast<std::size_t>(n_obs));
  for (std::int64_t r = 0; r < n_obs; ++r) all_rows[static_cast<std::size_t>(r)] = r;
  std::vector<double> column_mean(static_cast<std::size_t>(p));
  for (std::int64_t j = 0; j < p; ++j) {
    column_mean[static_cast<std::size_t>(j)] =
        r_mean(block.data() + j * n_obs, all_rows, false, nullptr);
  }
  for (std::int64_t j = 0; j < p; ++j) {
    const double* column_j = block.data() + j * n_obs;
    for (std::int64_t l = 0; l <= j; ++l) {
      const double* column_l = block.data() + l * n_obs;
      long double sum = 0.0L;
      for (std::int64_t r = 0; r < n_obs; ++r) {
        const double centred_j = column_j[r] - column_mean[static_cast<std::size_t>(j)];
        const double centred_l = column_l[r] - column_mean[static_cast<std::size_t>(l)];
        const double product = centred_j * centred_l;
        sum += product;
      }
      const double value = static_cast<double>(sum / static_cast<long double>(n_obs - 1));
      covariance[l + j * p] = value;
      covariance[j + l * p] = value;
    }
  }
  return Status::success();
}

Status subset_column(const CscView& matrix, Span<const int> rows, int column, Span<double> values) {
  if (values.size != rows.size) {
    return Status::failure(StatusCode::invalid_argument, "values must have one entry per row");
  }
  return gather_column(matrix, rows, column, values);
}

Status column_nonzero_rows(const CscView& matrix, int column, Span<int> rows, std::int64_t* n_found) {
  if (column < 1 || column > matrix.n_cols) {
    return Status::failure(StatusCode::invalid_argument, "column index out of range");
  }
  const int j = column - 1;
  std::int64_t found = 0;
  for (int k = matrix.column_pointer[j]; k < matrix.column_pointer[j + 1]; ++k) {
    // R's which(x != 0) drops a missing value rather than keeping it
    if (matrix.values[k] == 0.0 || std::isnan(matrix.values[k])) continue;
    if (rows.size != 0) {
      if (found >= rows.size) return Status::failure(StatusCode::invalid_argument, "rows buffer too small");
      rows[found] = matrix.row_index[k] + 1;
    }
    found += 1;
  }
  *n_found = found;
  return Status::success();
}

Status pair_variance_pratt(Span<const double> sigma, Span<const double> u, std::int64_t n_pairs,
                           std::int64_t n_genes, Span<double> v_pair, Span<double> v_diag,
                           Span<double> share_pct, Span<double> diag_pct, Span<double> sign,
                           double* v_total, double* v_diag_total, double* cross_cov_pct,
                           int* n_negative) {
  if (sigma.size != n_pairs * n_pairs || u.size != n_pairs * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "sigma must be n_pairs x n_pairs and u n_pairs x n_genes");
  }
  if (v_pair.size != n_pairs || v_diag.size != n_pairs || share_pct.size != n_pairs ||
      diag_pct.size != n_pairs || sign.size != n_pairs) {
    return Status::failure(StatusCode::invalid_argument, "outputs must have n_pairs entries");
  }
  std::vector<long double> pair_sum(static_cast<std::size_t>(n_pairs), 0.0L);
  std::vector<long double> square_sum(static_cast<std::size_t>(n_pairs), 0.0L);
  std::vector<double> slope(static_cast<std::size_t>(n_pairs));
  for (std::int64_t gene = 0; gene < n_genes; ++gene) {
    for (std::int64_t k = 0; k < n_pairs; ++k) {
      const double value = u[k + gene * n_pairs];
      slope[static_cast<std::size_t>(k)] = std::isfinite(value) ? value : 0.0;
    }
    for (std::int64_t k = 0; k < n_pairs; ++k) {
      double weighted = 0.0;
      for (std::int64_t l = 0; l < n_pairs; ++l) {
        weighted += sigma[k + l * n_pairs] * slope[static_cast<std::size_t>(l)];
      }
      pair_sum[static_cast<std::size_t>(k)] += slope[static_cast<std::size_t>(k)] * weighted;
      square_sum[static_cast<std::size_t>(k)] +=
          slope[static_cast<std::size_t>(k)] * slope[static_cast<std::size_t>(k)];
    }
  }
  long double total = 0.0L;
  long double diag_total = 0.0L;
  int negative = 0;
  bool any_missing = false;
  for (std::int64_t k = 0; k < n_pairs; ++k) {
    v_pair[k] = static_cast<double>(pair_sum[static_cast<std::size_t>(k)]);
    v_diag[k] = static_cast<double>(square_sum[static_cast<std::size_t>(k)]) * sigma[k + k * n_pairs];
    total += v_pair[k];
    diag_total += v_diag[k];
    // R's sum(V_pair < 0) is NA as soon as one comparison is NA
    if (std::isnan(v_pair[k])) any_missing = true;
    if (v_pair[k] < 0) negative += 1;
  }
  const double pratt_total = static_cast<double>(total);
  const double diagonal_total = static_cast<double>(diag_total);
  for (std::int64_t k = 0; k < n_pairs; ++k) {
    share_pct[k] = 100 * v_pair[k] / pratt_total;
    diag_pct[k] = 100 * v_diag[k] / diagonal_total;
    sign[k] = v_pair[k] > 0 ? 1.0 : (v_pair[k] < 0 ? -1.0 : (std::isnan(v_pair[k]) ? kQuietNaN : 0.0));
  }
  *v_total = pratt_total;
  *v_diag_total = diagonal_total;
  *cross_cov_pct = 100 * (pratt_total - diagonal_total) / pratt_total;
  *n_negative = any_missing ? -1 : negative;
  return Status::success();
}

}  // namespace pace
