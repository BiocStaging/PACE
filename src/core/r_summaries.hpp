// r_summaries.hpp -- R's own summary arithmetic, reproduced term for term so the
// compiled core returns what the R implementation returned.
//
// R accumulates sums in long double and rounds once at the end (src/main/summary.c);
// its mean adds a refinement pass (src/library/stats/src/cov.c uses the same one),
// and stats::var divides the centred sum of squares by n - 1. `skip_missing`
// reproduces na.rm = TRUE, which drops NA and NaN (both are NaN in C++) before
// anything is computed.
#ifndef PACE_R_SUMMARIES_HPP
#define PACE_R_SUMMARIES_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace pace {

// R's pmax(x, floor) for one value: NA and NaN stay missing.
inline double r_pmax(double value, double floor_value) {
  if (std::isnan(value)) return value;
  return value > floor_value ? value : floor_value;
}

// R's sum(x, na.rm = skip_missing) over the values at `rows`.
inline long double r_sum(const double* values, const std::vector<std::int64_t>& rows,
                         bool skip_missing) {
  long double sum = 0.0L;
  for (std::int64_t row : rows) {
    const double value = values[row];
    if (skip_missing && std::isnan(value)) continue;
    sum += value;
  }
  return sum;
}

// R's mean(x, na.rm = skip_missing) over the values at `rows`: the long double
// sum divided by the count, then refined by the long double sum of deviations.
// `n_used` receives the number of values that entered the mean.
inline double r_mean(const double* values, const std::vector<std::int64_t>& rows, bool skip_missing,
                     std::int64_t* n_used) {
  long double sum = 0.0L;
  std::int64_t used = 0;
  for (std::int64_t row : rows) {
    const double value = values[row];
    if (skip_missing && std::isnan(value)) continue;
    sum += value;
    used += 1;
  }
  if (n_used != nullptr) *n_used = used;
  if (used == 0) return std::numeric_limits<double>::quiet_NaN();
  long double mean = sum / static_cast<long double>(used);
  if (std::isfinite(static_cast<double>(mean))) {
    long double correction = 0.0L;
    for (std::int64_t row : rows) {
      const double value = values[row];
      if (skip_missing && std::isnan(value)) continue;
      correction += (value - mean);
    }
    mean = mean + correction / static_cast<long double>(used);
  }
  return static_cast<double>(mean);
}

// R's stats::var(x, na.rm = skip_missing) given the mean and count r_mean returned.
// Fewer than two values give NaN, which the binding reports as R's NA.
inline double r_variance(const double* values, const std::vector<std::int64_t>& rows, double mean,
                         std::int64_t n_used, bool skip_missing) {
  if (n_used < 2) return std::numeric_limits<double>::quiet_NaN();
  long double sum = 0.0L;
  for (std::int64_t row : rows) {
    const double value = values[row];
    if (skip_missing && std::isnan(value)) continue;
    const double centred = value - mean;
    const double square = centred * centred;
    sum += square;
  }
  return static_cast<double>(sum / static_cast<long double>(n_used - 1));
}

// R's mean(x, na.rm = skip_missing) over a contiguous array.
inline double r_mean_array(const double* values, std::int64_t n, bool skip_missing,
                           std::int64_t* n_used) {
  long double sum = 0.0L;
  std::int64_t used = 0;
  for (std::int64_t i = 0; i < n; ++i) {
    if (skip_missing && std::isnan(values[i])) continue;
    sum += values[i];
    used += 1;
  }
  if (n_used != nullptr) *n_used = used;
  if (used == 0) return std::numeric_limits<double>::quiet_NaN();
  long double mean = sum / static_cast<long double>(used);
  if (std::isfinite(static_cast<double>(mean))) {
    long double correction = 0.0L;
    for (std::int64_t i = 0; i < n; ++i) {
      if (skip_missing && std::isnan(values[i])) continue;
      correction += (values[i] - mean);
    }
    mean = mean + correction / static_cast<long double>(used);
  }
  return static_cast<double>(mean);
}

// R's stats::var(x, na.rm = skip_missing) over a contiguous array.
inline double r_variance_array(const double* values, std::int64_t n, bool skip_missing) {
  std::int64_t used = 0;
  const double mean = r_mean_array(values, n, skip_missing, &used);
  if (used < 2) return std::numeric_limits<double>::quiet_NaN();
  long double sum = 0.0L;
  for (std::int64_t i = 0; i < n; ++i) {
    if (skip_missing && std::isnan(values[i])) continue;
    const double centred = values[i] - mean;
    const double square = centred * centred;
    sum += square;
  }
  return static_cast<double>(sum / static_cast<long double>(used - 1));
}

// R's median(x, na.rm = TRUE): the middle value, or the mean of the two middle
// values for an even count. `scratch` is reordered in place.
inline double r_median(std::vector<double>& scratch) {
  std::vector<double>::iterator last = scratch.end();
  for (std::vector<double>::iterator it = scratch.begin(); it != last;) {
    if (std::isnan(*it)) {
      --last;
      std::iter_swap(it, last);
    } else {
      ++it;
    }
  }
  const std::int64_t n = static_cast<std::int64_t>(last - scratch.begin());
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  const std::int64_t half = (n + 1) / 2;
  std::nth_element(scratch.begin(), scratch.begin() + (half - 1), last);
  const double lower = scratch[static_cast<std::size_t>(half - 1)];
  if (n % 2 == 1) return lower;
  const double upper = *std::min_element(scratch.begin() + half, last);
  return static_cast<double>((static_cast<long double>(lower) + upper) / 2);
}

}  // namespace pace

#endif  // PACE_R_SUMMARIES_HPP
