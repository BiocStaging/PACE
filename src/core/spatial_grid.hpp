// spatial_grid.hpp -- fixed-radius neighbour search on a uniform grid.
//
// Reproduces dbscan::frNN(x, eps) semantics for the cells of one group:
//   * neighbours of cell i are all cells j != i of the same group with
//     dx*dx + dy*dy <= eps*eps (d == eps included; duplicates at d = 0 included);
//   * the distance is d = sqrt(dx*dx + dy*dy) with dx = x_i - x_j, dy = y_i - y_j,
//     each operation in its own statement (no FMA contraction).
// No neighbour list is stored: the caller supplies a visitor.
#ifndef PACE_SPATIAL_GRID_HPP
#define PACE_SPATIAL_GRID_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "core_types.hpp"

namespace pace {

class SpatialGrid {
 public:
  // x, y: coordinates of all n cells (caller-owned, alive while the grid is used).
  // cells: global indices (0-based) of the cells in this group.
  // radius: search radius, finite and > 0, in the units of x and y (callers
  // validate it; a degenerate radius yields a single bin and no neighbours).
  SpatialGrid(Span<const double> x, Span<const double> y,
              std::vector<std::int32_t> cells, double radius)
      : x_(x), y_(y), cells_(std::move(cells)), radius_(radius) {
    radius_sq_ = radius_ * radius_;
    // Query ranges use a radius widened by a few ulps so a pair whose rounded
    // squared distance equals radius^2 can never fall outside the scanned bins;
    // the exact test d2 <= radius^2 below still decides membership.
    query_radius_ = radius_ * (1.0 + 16.0 * std::numeric_limits<double>::epsilon());
    build();
  }

  const std::vector<std::int32_t>& cells() const { return cells_; }

  // Call visit(j, d) for every neighbour j of cell i (i must belong to this group).
  template <typename Visitor>
  void for_each_neighbour(std::int32_t i, Visitor&& visit) const {
    if (cells_.empty()) return;
    const double xi = x_[i];
    const double yi = y_[i];
    const std::int64_t bx_lo = bin_index(xi - query_radius_, x_min_, n_bins_x_);
    const std::int64_t bx_hi = bin_index(xi + query_radius_, x_min_, n_bins_x_);
    const std::int64_t by_lo = bin_index(yi - query_radius_, y_min_, n_bins_y_);
    const std::int64_t by_hi = bin_index(yi + query_radius_, y_min_, n_bins_y_);
    for (std::int64_t by = by_lo; by <= by_hi; ++by) {
      for (std::int64_t bx = bx_lo; bx <= bx_hi; ++bx) {
        const std::int64_t bin = by * n_bins_x_ + bx;
        for (std::int64_t k = bin_start_[bin]; k < bin_start_[bin + 1]; ++k) {
          const std::int32_t j = bin_cells_[k];
          if (j == i) continue;
          const double dx = xi - x_[j];
          const double dy = yi - y_[j];
          const double dx_sq = dx * dx;
          const double dy_sq = dy * dy;
          const double d_sq = dx_sq + dy_sq;
          if (d_sq <= radius_sq_) {
            const double d = std::sqrt(d_sq);
            visit(j, d);
          }
        }
      }
    }
  }

 private:
  // Bin of coordinate value v along an axis starting at `origin` with bin side
  // `bin_side_`, clamped to [0, n_bins - 1].
  std::int64_t bin_index(double v, double origin, std::int64_t n_bins) const {
    if (n_bins <= 1) return 0;
    const double position = std::floor((v - origin) / bin_side_);
    if (!(position > 0.0)) return 0;  // also catches NaN
    if (position >= static_cast<double>(n_bins - 1)) return n_bins - 1;
    return static_cast<std::int64_t>(position);
  }

  // Counting sort of the group's cells into bins. Bin side is the radius, widened
  // only when the extent would otherwise need more than ~2 bins per cell.
  void build() {
    if (cells_.empty()) return;
    x_min_ = x_[cells_[0]];
    double x_max = x_min_;
    y_min_ = y_[cells_[0]];
    double y_max = y_min_;
    for (std::int32_t c : cells_) {
      x_min_ = std::min(x_min_, x_[c]);
      x_max = std::max(x_max, x_[c]);
      y_min_ = std::min(y_min_, y_[c]);
      y_max = std::max(y_max, y_[c]);
    }
    const double width = x_max - x_min_;
    const double height = y_max - y_min_;
    const double max_bins = 2.0 * static_cast<double>(cells_.size()) + 16.0;
    bin_side_ = radius_;
    // Doubling keeps the bin count bounded for any extent, including very long,
    // thin sections; a larger bin side only adds candidates, never misses pairs.
    // The loop is capped (2^1100 exceeds any finite double) so it cannot spin.
    const bool usable_geometry = std::isfinite(bin_side_) && bin_side_ > 0.0 &&
                                 std::isfinite(width) && std::isfinite(height);
    if (usable_geometry) {
      for (int doubling = 0; doubling < 1100; ++doubling) {
        const double bins = (std::floor(width / bin_side_) + 1.0) * (std::floor(height / bin_side_) + 1.0);
        if (bins <= max_bins) break;
        bin_side_ *= 2.0;
      }
    }
    if (!usable_geometry || !std::isfinite(bin_side_)) {
      // Degenerate input (callers validate radius and coordinates first): a single bin.
      bin_side_ = std::numeric_limits<double>::infinity();
      n_bins_x_ = 1;
      n_bins_y_ = 1;
    } else {
      n_bins_x_ = static_cast<std::int64_t>(std::floor(width / bin_side_)) + 1;
      n_bins_y_ = static_cast<std::int64_t>(std::floor(height / bin_side_)) + 1;
    }

    const std::int64_t n_bins = n_bins_x_ * n_bins_y_;
    bin_start_.assign(n_bins + 1, 0);
    std::vector<std::int64_t> cell_bin(cells_.size());
    for (std::size_t k = 0; k < cells_.size(); ++k) {
      const std::int32_t c = cells_[k];
      const std::int64_t bx = bin_index(x_[c], x_min_, n_bins_x_);
      const std::int64_t by = bin_index(y_[c], y_min_, n_bins_y_);
      cell_bin[k] = by * n_bins_x_ + bx;
      bin_start_[cell_bin[k] + 1] += 1;
    }
    for (std::int64_t b = 0; b < n_bins; ++b) bin_start_[b + 1] += bin_start_[b];
    bin_cells_.assign(cells_.size(), 0);
    std::vector<std::int64_t> fill_position(bin_start_.begin(), bin_start_.end() - 1);
    for (std::size_t k = 0; k < cells_.size(); ++k) {
      bin_cells_[fill_position[cell_bin[k]]++] = cells_[k];
    }
  }

  Span<const double> x_;
  Span<const double> y_;
  std::vector<std::int32_t> cells_;
  double radius_ = 0.0;
  double radius_sq_ = 0.0;
  double query_radius_ = 0.0;
  double x_min_ = 0.0;
  double y_min_ = 0.0;
  double bin_side_ = 1.0;
  std::int64_t n_bins_x_ = 1;
  std::int64_t n_bins_y_ = 1;
  std::vector<std::int64_t> bin_start_;
  std::vector<std::int32_t> bin_cells_;
};

}  // namespace pace

#endif  // PACE_SPATIAL_GRID_HPP
