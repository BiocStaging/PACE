// neighbourhood.cpp -- implementation of neighbourhood.hpp.
#include "fp_no_contract.hpp"  // must precede the distance and kernel arithmetic

#include "neighbourhood.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

#include "spatial_grid.hpp"
#include "thread_pool.hpp"

namespace pace {

namespace {

constexpr double kPi = 3.141592653589793238462643383280;  // R's `pi` (M_PI)
constexpr std::int64_t kCellBlock = 256;                   // cells per parallel block

// Grids for the searchable groups and, per cell, the grid it belongs to (-1: none).
struct GroupGrids {
  std::vector<SpatialGrid> grids;
  std::vector<int> grid_of_cell;
};

// Global indices of the cells in each group code 0..n_groups-1, ascending.
std::vector<std::vector<std::int32_t>> cells_by_group(Span<const int> group, int n_groups) {
  std::vector<std::vector<std::int32_t>> members(n_groups);
  for (std::int64_t i = 0; i < group.size; ++i) {
    const int g = group[i];
    if (g >= 0 && g < n_groups) members[g].push_back(static_cast<std::int32_t>(i));
  }
  return members;
}

// Largest group code + 1 (0 when every code is negative).
int count_groups(Span<const int> group) {
  int n_groups = 0;
  for (std::int64_t i = 0; i < group.size; ++i) n_groups = std::max(n_groups, group[i] + 1);
  return n_groups;
}

// Build search grids. per_group = false: one grid over all cells. per_group = true:
// one grid per group with at least `min_cells` cells; smaller groups get none.
GroupGrids build_grids(Span<const double> x, Span<const double> y, Span<const int> group,
                       bool per_group, double radius, std::int64_t min_cells) {
  GroupGrids out;
  const std::int64_t n = x.size;
  out.grid_of_cell.assign(n, -1);
  if (!per_group) {
    std::vector<std::int32_t> all_cells(n);
    for (std::int64_t i = 0; i < n; ++i) all_cells[i] = static_cast<std::int32_t>(i);
    if (n >= min_cells) {
      out.grids.emplace_back(x, y, std::move(all_cells), radius);
      std::fill(out.grid_of_cell.begin(), out.grid_of_cell.end(), 0);
    }
    return out;
  }
  const int n_groups = count_groups(group);
  std::vector<std::vector<std::int32_t>> members = cells_by_group(group, n_groups);
  for (int g = 0; g < n_groups; ++g) {
    if (static_cast<std::int64_t>(members[g].size()) < min_cells) continue;
    const int grid_id = static_cast<int>(out.grids.size());
    for (std::int32_t c : members[g]) out.grid_of_cell[c] = grid_id;
    out.grids.emplace_back(x, y, std::move(members[g]), radius);
  }
  return out;
}

Status check_same_length(std::int64_t expected, std::int64_t actual, const char* name) {
  if (expected == actual) return Status::success();
  return Status::failure(StatusCode::invalid_argument,
                         std::string("length mismatch for ") + name);
}

// Distance from (x0, y0) to the rectangle edge along one angle, as R's
// pmin(ifelse(cos > 0, (xmax - x0)/cos, Inf), ifelse(cos < 0, (xmin - x0)/cos, Inf),
//      ifelse(sin > 0, (ymax - y0)/sin, Inf), ifelse(sin < 0, (ymin - y0)/sin, Inf)).
double edge_distance(double x0, double y0, double cos_t, double sin_t,
                     double x_min, double x_max, double y_min, double y_max) {
  const double inf = std::numeric_limits<double>::infinity();
  const double d_right = cos_t > 0.0 ? (x_max - x0) / cos_t : inf;
  const double d_left = cos_t < 0.0 ? (x_min - x0) / cos_t : inf;
  const double d_top = sin_t > 0.0 ? (y_max - y0) / sin_t : inf;
  const double d_bottom = sin_t < 0.0 ? (y_min - y0) / sin_t : inf;
  double d_max = d_right;
  if (d_left < d_max) d_max = d_left;
  if (d_top < d_max) d_max = d_top;
  if (d_bottom < d_max) d_max = d_bottom;
  return d_max;
}

// Quadrature fraction given sum_theta r_eff^2 (long double, as R's sum()).
double fraction_from_sum(long double sum_sq, double r, std::int64_t n_angles) {
  const double sum_double = static_cast<double>(sum_sq);
  const double numerator = (kPi * sum_double) / static_cast<double>(n_angles);
  const double norm_factor = kPi * (r * r);
  return numerator / norm_factor;
}

// Edge fraction of one cell by full quadrature.
double cell_area_fraction(double x0, double y0, double r, double x_min, double x_max,
                          double y_min, double y_max, Span<const double> cos_theta,
                          Span<const double> sin_theta) {
  long double sum_sq = 0.0L;
  for (std::int64_t k = 0; k < cos_theta.size; ++k) {
    const double d_max = edge_distance(x0, y0, cos_theta[k], sin_theta[k], x_min, x_max, y_min, y_max);
    const double r_eff = d_max < r ? d_max : r;  // R's pmin(r, d_max)
    const double r_eff_sq = r_eff * r_eff;
    sum_sq += r_eff_sq;
  }
  return fraction_from_sum(sum_sq, r, cos_theta.size);
}

// Edge fraction shared by every cell at least r from all four edges.
double interior_area_fraction(double r, std::int64_t n_angles) {
  long double sum_sq = 0.0L;
  const double r_sq = r * r;
  for (std::int64_t k = 0; k < n_angles; ++k) sum_sq += r_sq;
  return fraction_from_sum(sum_sq, r, n_angles);
}

// True when every ray from (x0, y0) reaches the rectangle edge at distance >= r,
// so r_eff = r at every angle. (x_max - x0) / cos >= x_max - x0 for 0 < cos <= 1,
// and the same holds for the other three edges.
bool is_interior(double x0, double y0, double r, double x_min, double x_max,
                 double y_min, double y_max) {
  return (x_max - x0) >= r && -(x_min - x0) >= r && (y_max - y0) >= r && -(y_min - y0) >= r;
}

}  // namespace

Status neighbour_kernels(Span<const double> x, Span<const double> y,
                         Span<const int> neighbour_type, int n_types,
                         Span<const int> group, bool per_group,
                         double h_bio, double h_tech, double eps,
                         Span<double> k_bio, Span<double> k_tech,
                         int n_threads, const InterruptCheck& interrupted) {
  const std::int64_t n = x.size;
  Status status = check_same_length(n, y.size, "y");
  if (!status.is_ok()) return status;
  status = check_same_length(n, neighbour_type.size, "neighbour_type");
  if (!status.is_ok()) return status;
  status = check_same_length(n, group.size, "group");
  if (!status.is_ok()) return status;
  status = check_same_length(n * n_types, k_bio.size, "k_bio");
  if (!status.is_ok()) return status;
  status = check_same_length(n * n_types, k_tech.size, "k_tech");
  if (!status.is_ok()) return status;
  if (!(eps > 0.0)) return Status::failure(StatusCode::invalid_argument, "eps must be positive");

  const GroupGrids grids = build_grids(x, y, group, per_group, eps, 2);
  const double h_bio_sq = h_bio * h_bio;

  struct Neighbour {
    double distance;
    int type;
  };

  auto body = [&](std::int64_t begin, std::int64_t end) {
    std::vector<Neighbour> neighbours;
    for (std::int64_t i = begin; i < end; ++i) {
      const int grid_id = grids.grid_of_cell[i];
      if (grid_id < 0) continue;
      neighbours.clear();
      grids.grids[grid_id].for_each_neighbour(static_cast<std::int32_t>(i), [&](std::int32_t j, double d) {
        const int type = neighbour_type[j];
        if (type >= 0 && type < n_types) neighbours.push_back(Neighbour{d, type});
      });
      // Increasing distance: the frNN order. Equal distances give equal kernel
      // values, so their relative order cannot change a sum.
      std::sort(neighbours.begin(), neighbours.end(),
                [](const Neighbour& a, const Neighbour& b) { return a.distance < b.distance; });
      for (const Neighbour& nb : neighbours) {
        const double tech_weight = std::exp(-nb.distance / h_tech);
        const double distance_sq = nb.distance * nb.distance;
        const double bio_weight = std::exp(-distance_sq / h_bio_sq);
        const std::int64_t slot = i + static_cast<std::int64_t>(nb.type) * n;
        k_tech[slot] += tech_weight;
        k_bio[slot] += bio_weight;
      }
    }
  };
  return parallel_for(n, n_threads, kCellBlock, body, interrupted);
}

Status neighbour_counts(Span<const double> x, Span<const double> y,
                        Span<const int> group, bool per_group, double eps,
                        Span<int> counts, int n_threads, const InterruptCheck& interrupted) {
  const std::int64_t n = x.size;
  Status status = check_same_length(n, counts.size, "counts");
  if (!status.is_ok()) return status;
  status = check_same_length(n, group.size, "group");
  if (!status.is_ok()) return status;
  const GroupGrids grids = build_grids(x, y, group, per_group, eps, 2);
  auto body = [&](std::int64_t begin, std::int64_t end) {
    for (std::int64_t i = begin; i < end; ++i) {
      int count = 0;
      const int grid_id = grids.grid_of_cell[i];
      if (grid_id >= 0) {
        grids.grids[grid_id].for_each_neighbour(static_cast<std::int32_t>(i),
                                                [&](std::int32_t, double) { count += 1; });
      }
      counts[i] = count;
    }
  };
  return parallel_for(n, n_threads, kCellBlock, body, interrupted);
}

Status area_fractions(Span<const double> x, Span<const double> y, Span<const int> rows,
                      double r, double x_min, double x_max, double y_min, double y_max,
                      Span<const double> cos_theta, Span<const double> sin_theta,
                      Span<double> fraction, int n_threads, const InterruptCheck& interrupted) {
  Status status = check_same_length(rows.size, fraction.size, "fraction");
  if (!status.is_ok()) return status;
  status = check_same_length(cos_theta.size, sin_theta.size, "sin_theta");
  if (!status.is_ok()) return status;
  if (cos_theta.size < 1) return Status::failure(StatusCode::invalid_argument, "no angles");
  const double interior_value = interior_area_fraction(r, cos_theta.size);
  auto body = [&](std::int64_t begin, std::int64_t end) {
    for (std::int64_t k = begin; k < end; ++k) {
      const int cell = rows[k];
      const double x0 = x[cell];
      const double y0 = y[cell];
      if (is_interior(x0, y0, r, x_min, x_max, y_min, y_max)) {
        fraction[k] = interior_value;
      } else {
        fraction[k] = cell_area_fraction(x0, y0, r, x_min, x_max, y_min, y_max, cos_theta, sin_theta);
      }
    }
  };
  return parallel_for(rows.size, n_threads, kCellBlock, body, interrupted);
}

// ---------------------------------------------------------------------------
// AmbientFieldBuilder
// ---------------------------------------------------------------------------

struct AmbientFieldBuilder::State {
  Span<const double> x;
  Span<const double> y;
  Span<const int> type_code;
  double radius = 0.0;
  double h_tech = 0.0;
  GroupGrids grids;
  std::vector<std::int64_t> column_pointer;  // length n + 1
};

AmbientFieldBuilder::AmbientFieldBuilder() : state_(new State()) {}
AmbientFieldBuilder::~AmbientFieldBuilder() = default;

std::int64_t AmbientFieldBuilder::nnz() const {
  return state_->column_pointer.empty() ? 0 : state_->column_pointer.back();
}

const std::vector<std::int64_t>& AmbientFieldBuilder::column_pointer() const {
  return state_->column_pointer;
}

Status AmbientFieldBuilder::prepare(Span<const double> x, Span<const double> y,
                                    Span<const int> type_code, Span<const int> image,
                                    int n_images, double h_tech, bool edge_correct,
                                    Span<const double> cos_theta, Span<const double> sin_theta,
                                    Span<double> edge_fraction, int n_threads,
                                    const InterruptCheck& interrupted) {
  const std::int64_t n = x.size;
  Status status = check_same_length(n, y.size, "y");
  if (!status.is_ok()) return status;
  status = check_same_length(n, type_code.size, "type_code");
  if (!status.is_ok()) return status;
  status = check_same_length(n, image.size, "image");
  if (!status.is_ok()) return status;
  status = check_same_length(n, edge_fraction.size, "edge_fraction");
  if (!status.is_ok()) return status;

  State& s = *state_;
  s.x = x;
  s.y = y;
  s.type_code = type_code;
  s.h_tech = h_tech;
  s.radius = 3.0 * h_tech;

  // Edge fractions against each image's bounding box, one parallel pass over cells.
  std::fill(edge_fraction.data, edge_fraction.data + n, 1.0);
  if (edge_correct) {
    status = check_same_length(cos_theta.size, sin_theta.size, "sin_theta");
    if (!status.is_ok()) return status;
    if (cos_theta.size < 1) return Status::failure(StatusCode::invalid_argument, "no angles");
    const std::vector<std::vector<std::int32_t>> members = cells_by_group(image, n_images);
    std::vector<double> x_min(n_images), x_max(n_images), y_min(n_images), y_max(n_images);
    for (int g = 0; g < n_images; ++g) {
      const std::vector<std::int32_t>& rows = members[g];
      if (rows.empty()) continue;
      x_min[g] = x[rows[0]];
      x_max[g] = x_min[g];
      y_min[g] = y[rows[0]];
      y_max[g] = y_min[g];
      for (std::int32_t c : rows) {
        x_min[g] = std::min(x_min[g], x[c]);
        x_max[g] = std::max(x_max[g], x[c]);
        y_min[g] = std::min(y_min[g], y[c]);
        y_max[g] = std::max(y_max[g], y[c]);
      }
    }
    const double interior_value = interior_area_fraction(s.radius, cos_theta.size);
    auto edge_body = [&](std::int64_t begin, std::int64_t end) {
      for (std::int64_t i = begin; i < end; ++i) {
        const int g = image[i];
        if (g < 0 || g >= n_images) continue;
        if (is_interior(x[i], y[i], s.radius, x_min[g], x_max[g], y_min[g], y_max[g])) {
          edge_fraction[i] = interior_value;
        } else {
          edge_fraction[i] = cell_area_fraction(x[i], y[i], s.radius, x_min[g], x_max[g],
                                                y_min[g], y_max[g], cos_theta, sin_theta);
        }
      }
    };
    status = parallel_for(n, n_threads, kCellBlock, edge_body, interrupted);
    if (!status.is_ok()) return status;
  }

  // Column counts: heterotypic neighbours of each cell within its image.
  s.grids = build_grids(x, y, image, true, s.radius, 2);
  std::vector<std::int64_t> counts(n, 0);
  auto count_body = [&](std::int64_t begin, std::int64_t end) {
    for (std::int64_t j = begin; j < end; ++j) {
      const int grid_id = s.grids.grid_of_cell[j];
      if (grid_id < 0) continue;
      std::int64_t count = 0;
      const int type_j = type_code[j];
      s.grids.grids[grid_id].for_each_neighbour(static_cast<std::int32_t>(j), [&](std::int32_t i, double) {
        if (type_code[i] != type_j) count += 1;
      });
      counts[j] = count;
    }
  };
  status = parallel_for(n, n_threads, kCellBlock, count_body, interrupted);
  if (!status.is_ok()) return status;

  s.column_pointer.assign(n + 1, 0);
  for (std::int64_t j = 0; j < n; ++j) s.column_pointer[j + 1] = s.column_pointer[j] + counts[j];
  if (nnz() > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
    return Status::failure(StatusCode::invalid_argument,
                           "ambient field has more non-zeros than a dgCMatrix can hold");
  }
  return Status::success();
}

Status AmbientFieldBuilder::fill(Span<int> row_index, Span<double> values,
                                 Span<const double> edge_fraction, int n_threads,
                                 const InterruptCheck& interrupted) {
  State& s = *state_;
  Status status = check_same_length(nnz(), row_index.size, "row_index");
  if (!status.is_ok()) return status;
  status = check_same_length(nnz(), values.size, "values");
  if (!status.is_ok()) return status;
  const std::int64_t n = s.x.size;

  struct Entry {
    std::int32_t row;
    double distance;
  };

  auto fill_body = [&](std::int64_t begin, std::int64_t end) {
    std::vector<Entry> entries;
    for (std::int64_t j = begin; j < end; ++j) {
      const int grid_id = s.grids.grid_of_cell[j];
      if (grid_id < 0) continue;
      entries.clear();
      const int type_j = s.type_code[j];
      s.grids.grids[grid_id].for_each_neighbour(static_cast<std::int32_t>(j), [&](std::int32_t i, double d) {
        if (s.type_code[i] != type_j) entries.push_back(Entry{i, d});
      });
      std::sort(entries.begin(), entries.end(),
                [](const Entry& a, const Entry& b) { return a.row < b.row; });
      std::int64_t position = s.column_pointer[j];
      for (const Entry& entry : entries) {
        // d_ij is symmetric (dx and -dx square identically), so the value equals
        // R's exp(-d / h_tech) / af[i] computed from cell i's neighbour list.
        const double kernel_weight = std::exp(-entry.distance / s.h_tech);
        row_index[position] = entry.row;
        values[position] = kernel_weight / edge_fraction[entry.row];
        position += 1;
      }
    }
  };
  return parallel_for(n, n_threads, kCellBlock, fill_body, interrupted);
}

Status same_type_fraction(Span<const double> x, Span<const double> y,
                          Span<const int> type_code, Span<const int> image, int n_images,
                          double radius, int min_image_cells, Span<double> fraction,
                          int n_threads, const InterruptCheck& interrupted) {
  const std::int64_t n = x.size;
  Status status = check_same_length(n, fraction.size, "fraction");
  if (!status.is_ok()) return status;
  status = check_same_length(n, type_code.size, "type_code");
  if (!status.is_ok()) return status;
  status = check_same_length(n, image.size, "image");
  if (!status.is_ok()) return status;
  (void)n_images;
  const GroupGrids grids = build_grids(x, y, image, true, radius, min_image_cells);
  auto body = [&](std::int64_t begin, std::int64_t end) {
    for (std::int64_t i = begin; i < end; ++i) {
      const int grid_id = grids.grid_of_cell[i];
      if (grid_id < 0) {
        fraction[i] = std::numeric_limits<double>::quiet_NaN();
        continue;
      }
      std::int64_t n_neighbours = 0;
      std::int64_t n_same = 0;
      const int type_i = type_code[i];
      grids.grids[grid_id].for_each_neighbour(static_cast<std::int32_t>(i), [&](std::int32_t j, double) {
        n_neighbours += 1;
        if (type_code[j] == type_i) n_same += 1;
      });
      if (n_neighbours == 0) {
        fraction[i] = 1.0;
      } else {
        const long double ratio = static_cast<long double>(n_same) / static_cast<long double>(n_neighbours);
        fraction[i] = static_cast<double>(ratio);
      }
    }
  };
  return parallel_for(n, n_threads, kCellBlock, body, interrupted);
}

}  // namespace pace
