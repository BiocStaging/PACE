// dispersion.hpp -- the numerical optimisation the fit does outside the solve:
// the per-gene negative binomial dispersion MLE, and the prior degrees of
// freedom of the adaptive tau shrinkage.
//
// Both were R calls (`optimize()` over `dnbinom()`, and `uniroot()` over
// `trigamma()`). The algorithms are here -- Brent's function minimiser and
// Brent's root finder, in the form R's optimize() and uniroot() use, plus the
// likelihood assembly and the loop over genes. The two elementary special
// functions are NOT reimplemented: the caller passes them in, and the R binding
// passes R's own, so a fit gives the same numbers to the last bit while the core
// stays free of R. A caller outside R supplies its own pair.
#ifndef PACE_DISPERSION_HPP
#define PACE_DISPERSION_HPP

#include <cstdint>

#include "core_types.hpp"
#include "count_stats.hpp"

namespace pace {

// log f(x; size, mu) of the negative binomial in its mean parameterisation:
// R's dnbinom(x, size = size, mu = mu, log = TRUE).
using LogDensity = double (*)(double x, double size, double mu);

// The trigamma function, R's trigamma().
using Trigamma = double (*)(double x);

// R's optimize(): Brent's parabolic-interpolation minimiser on [lower, upper],
// stopping at eps |x| + tol/3 with eps the square root of the machine epsilon.
// Returns the minimiser; `n_evaluations` (may be null) counts the objective calls.
double brent_minimise(double lower, double upper, double tolerance,
                      const std::function<double(double)>& objective,
                      std::int64_t* n_evaluations);

// R's uniroot(): Brent's zeroin on a bracketing interval. Returns the root, or
// NaN when the ends do not bracket a sign change (R raises an error there, which
// its callers catch). `max_iterations` matches R's default of 1000.
double brent_root(double lower, double upper, double tolerance, int max_iterations,
                  const std::function<double(double)>& objective);

// The dispersion MLE of one gene, given its counts and fitted means:
//   alpha = argmin_a -sum_i log f(y_i; size, mu_i),  size = mu_i / a (NB1) or 1 / a (NB2),
// minimised over log a in [-6, 4]. Fewer than ten cells gives NaN. See
// dispersion_chunk() for `zero_collapse` and `max_cells`.
double dispersion_mle(Span<const double> counts, Span<const double> mu, bool nb2,
                      bool zero_collapse, double max_cells, LogDensity density);

// The dispersion MLE of every gene in one chunk.
//
// For gene j the fitted mean is rebuilt from the chunk's linear predictor and
// the ambient field, as the solver's own final pass does, and then
//
//   alpha_j = argmin_a  -sum_i log f(y_ij; size, mu_ij),
//     size = mu_ij / a   (NB1)   or   size = 1 / a   (NB2)
//
// minimised over log a in [-6, 4]. Cells with a non-finite or vanishing fitted
// mean are dropped, and a gene left with fewer than ten is reported NaN.
//
// `zero_collapse` (NB1 only) uses log f(0; mu/a, mu) = -(mu/a) log1p(a) to sum
// the zero-count cells in closed form instead of calling the density on each.
// It is exact algebra, but it reassociates the sum, so it changes the last
// digits of alpha and is off unless the caller asks for it.
//
// `max_cells` subsamples evenly when a gene has more cells than that, as the R
// implementation's `alpha_max_n` did; a non-finite value keeps every cell.
//
// A count that is not a whole number has no negative binomial density: those
// genes are reported NaN and counted in `n_noninteger`, which the caller warns
// about, rather than calling into the density (which would want to warn from a
// worker thread).
//
// Shapes: `eta` is n * n_genes column-major; `offset` and `rho` have n entries;
// `alpha` has n_genes. The counts and the ambient field are read one gene at a
// time from their sparse blocks.
Status dispersion_chunk(Span<const double> eta, const GeneBlock& counts, const GeneBlock& ambient,
                        Span<const double> offset, Span<const double> rho, bool nb2,
                        bool zero_collapse, double max_cells, LogDensity density, std::int64_t n,
                        std::int64_t n_genes, Span<double> alpha, std::int64_t* n_noninteger,
                        int n_threads, const InterruptCheck& interrupted);

// The prior degrees of freedom of the adaptive tau shrinkage, from the
// cross-gene variance of log s2 and how many genes entered it:
//   excess = variance - trigamma(1/2);  d0 = min(2 x, d0_max) with trigamma(x) = excess
// With fewer than five genes, a non-finite variance, or no excess over the
// sampling term, the cap d0_max is returned, as the R implementation did.
double estimate_d0(double log_variance, std::int64_t n_used, double d0_max, Trigamma trigamma);

}  // namespace pace

#endif  // PACE_DISPERSION_HPP
