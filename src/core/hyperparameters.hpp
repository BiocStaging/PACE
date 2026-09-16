// hyperparameters.hpp -- the variance-component and weight arithmetic of the
// solver: the EM update of tau, its three shrinkage modes, the cap, and the
// data-informed prior weights.
//
// What stays in R is the estimation of the prior degrees of freedom d0, which
// needs R's own trigamma() and uniroot(); reimplementing those here would move
// the shrinkage in the last bits for no gain. The core returns the cross-gene
// log-variance d0 is estimated from, and takes d0 back.
#ifndef PACE_HYPERPARAMETERS_HPP
#define PACE_HYPERPARAMETERS_HPP

#include <cstdint>

#include "core_types.hpp"

namespace pace {

// The EM update of the variance components, per random-effect column:
//   tau_k = max(mean_g(u_kg^2 + v_kg), 1e-6)      (na.rm, R's refined mean)
// Shapes: `u` and `re_var` are q * n_genes column-major; `tau` has q entries.
Status tau_em_update(Span<const double> u, Span<const double> re_var, std::int64_t q,
                     std::int64_t n_genes, Span<double> tau);

// Hierarchical shrinkage of one block's tau towards its group-size-weighted
// mean, per term (Hierarchical mode):
//   lambda   = lambda_factor * median(n_group)
//   weight_c = n_c / (n_c + lambda)
//   global_t = sum_c n_c tau_tc / sum_c n_c
//   out_tc   = max(weight_c tau_tc + (1 - weight_c) global_t, 1e-6)
// A group size of zero or less is read as one, as the R implementation did.
// Shapes: `tau` and `out` are n_terms * n_groups column-major (terms vary
// fastest is NOT assumed: the matrix is R's tau[term, group]); `n_group` has
// n_groups entries. With a mismatched `n_group` the input is returned unchanged.
Status tau_hierarchical(Span<const double> tau, Span<const int> n_group, std::int64_t n_terms,
                        std::int64_t n_groups, double lambda_factor, Span<double> out);

// The per-row summaries the adaptive empirical-Bayes mode needs, over
// vals_kg = max(s2_kg, 1e-9) * reml_factor:
//   panel_k        = mean(vals_k)            (na.rm, R's refined mean)
//   panel_median_k = median(vals_k)          (na.rm, R's median)
//   log_variance_k = var(log(vals_k))        over the finite entries only
//   n_log_finite_k = how many entries entered that variance
// The caller turns log_variance into d0 with R's trigamma and uniroot.
Status tau_eb_summaries(Span<const double> s2, std::int64_t q, std::int64_t n_genes,
                        double reml_factor, Span<double> panel, Span<double> panel_median,
                        Span<double> log_variance, Span<int> n_log_finite);

// The adaptive empirical-Bayes shrinkage itself, given d0:
//   tau_kg = max((d0_k panel_k + vals_kg) / (d0_k + 1), floor_k)
//   floor_k = max(panel_median_k / 100, 1e-4), or `tau_floor` when it is finite
// Shapes: `s2` and `out` are q * n_genes; the per-row spans have q entries.
Status tau_eb_apply(Span<const double> s2, std::int64_t q, std::int64_t n_genes,
                    double reml_factor, Span<const double> panel,
                    Span<const double> panel_median, Span<const double> d0, double tau_floor,
                    Span<double> out);

// The half-Cauchy mode: per row, an ECM over
//   tau_g = max((1 / a_g + s2_g / 2) / 2, floor)
//   a_g   = max((1 / tau_g + 1 / lambda2) / 2, 1e-9)
//   lambda2 = max(2 mean(1 / a_g), 1e-4)
// started from lambda2 = max(median(s2_g), 1e-3) and a_g = lambda2, with
// floor = max(mean(s2_g) / 100, tau_floor) and s2_g = max(s2, 1e-9).
// `lambda2_prev` (q) and `a_prev` (q * n_genes) are warm starts from the
// previous iteration; either may be empty, and a non-finite entry falls back to
// the cold start. `lambda2_out`, `a_out` and `panel_out` may be empty.
Status tau_half_cauchy(Span<const double> s2, std::int64_t q, std::int64_t n_genes, int n_em_iter,
                       double tau_floor, Span<const double> lambda2_prev,
                       Span<const double> a_prev, Span<double> out, Span<double> lambda2_out,
                       Span<double> a_out, Span<double> panel_out);

// Cap the variance components, reporting how many bound (the caller warns).
Status tau_clamp(Span<double> tau, std::int64_t size, double tau_max, std::int64_t* n_binding);

// The data-informed prior weights, given each random-effect row's focal cell
// type and its neighbour-kernel scale:
//   W_kg = detection_rate[focal_k, g] * scale_k          (rows with a focal)
//   W_kg = 1                                             (rows without one)
// then each row is divided by its own maximum (so its most informative gene has
// weight 1) and floored at 1e-8.
// Shapes: `detection_rate` is n_focals * n_genes column-major; `focal_of_row`
// and `scale_of_row` have q entries, the focal being 1-based with 0 meaning
// "leave this row at one"; `weights` is q * n_genes.
Status data_informed_weights(Span<const double> detection_rate, std::int64_t n_focals,
                             std::int64_t n_genes, Span<const int> focal_of_row,
                             Span<const double> scale_of_row, std::int64_t q,
                             Span<double> weights);

}  // namespace pace

#endif  // PACE_HYPERPARAMETERS_HPP
