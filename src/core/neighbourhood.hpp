// neighbourhood.hpp -- neighbour kernels, ambient contamination field, edge
// correction and same-type neighbour fractions, computed without storing pairs.
//
// All functions reproduce the R reference implementations in R/pace-core.R
// (f55d976) bit for bit on the same platform: same distance arithmetic as
// dbscan::frNN, same accumulation order as Matrix::sparseMatrix, and R's
// long double sums where R uses them. Coordinates are in the caller's units
// (micrometres in PACE). Cell indices are 0-based.
#ifndef PACE_NEIGHBOURHOOD_HPP
#define PACE_NEIGHBOURHOOD_HPP

#include <cstdint>
#include <memory>
#include <vector>

#include "core_types.hpp"

namespace pace {

// Every function refuses (Status invalid_argument) non-finite coordinates and
// radii or bandwidths that are not finite and positive.

// Neighbour kernels, accumulated per neighbour cell type.
//
//   K_tech[i, t] = sum over neighbours j of type t of exp(-d_ij / h_tech)
//   K_bio[i, t]  = sum over neighbours j of type t of exp(-d_ij^2 / h_bio^2)
//
// Neighbours are cells j != i within `eps` (d == eps included), taken over all
// cells when `per_group` is false, or only within cell i's group when true (a
// group with fewer than 2 cells, or group -1, gives no neighbours). Cells with
// neighbour_type -1 are never counted as neighbours but still get output rows.
// Each (i, t) sum runs over neighbours in increasing distance, the order in which
// frNN lists them and Matrix::sparseMatrix adds duplicate entries.
//
// Shapes: x, y, neighbour_type, group have length n. k_bio and k_tech have length
// n * n_types (column-major n x n_types), are caller-owned and must be zero on entry.
Status neighbour_kernels(Span<const double> x, Span<const double> y,
                         Span<const int> neighbour_type, int n_types,
                         Span<const int> group, bool per_group,
                         double h_bio, double h_tech, double eps,
                         Span<double> k_bio, Span<double> k_tech,
                         int n_threads, const InterruptCheck& interrupted);

// Number of neighbours of every cell under the same rules as neighbour_kernels(),
// without type filtering (equals lengths(dbscan::frNN(...)$id)). Used by tests.
Status neighbour_counts(Span<const double> x, Span<const double> y,
                        Span<const int> group, bool per_group, double eps,
                        Span<int> counts, int n_threads, const InterruptCheck& interrupted);

// Neighbour lists under the same rules (test helper): for each cell i, in cell
// order, its neighbours j (0-based) sorted by index, with their distances.
// `offsets` gets n + 1 entries; `neighbours` and `distances` are appended.
Status neighbour_lists(Span<const double> x, Span<const double> y,
                       Span<const int> group, bool per_group, double eps,
                       std::vector<std::int64_t>& offsets, std::vector<int>& neighbours,
                       std::vector<double>& distances, int n_threads,
                       const InterruptCheck& interrupted);

// Isotropic edge correction for the cells `rows` against the rectangle
// [x_min, x_max] x [y_min, y_max]: the fraction of the disc of radius r around
// the cell that lies inside the rectangle, by quadrature over the angles whose
// cosines and sines are supplied (computed in R so the angles are identical):
//
//   r_eff(theta) = min(r, distance from the cell to the rectangle edge along theta)
//   fraction     = (pi * sum_theta r_eff^2 / n_angles) / (pi * r^2)
//
// A cell at least r from every edge has r_eff = r for every angle, so all such
// cells share one value, computed once. Output has length rows.size.
Status area_fractions(Span<const double> x, Span<const double> y, Span<const int> rows,
                      double r, double x_min, double x_max, double y_min, double y_max,
                      Span<const double> cos_theta, Span<const double> sin_theta,
                      Span<double> fraction, int n_threads, const InterruptCheck& interrupted);

// Sparse cross-cell-type ambient weight matrix W (n x n, compressed column form):
//
//   W[i, j] = exp(-d_ij / h_tech) / edge_fraction[i]
//
// for cells i, j of the same image, of different type codes, within 3 * h_tech.
// edge_fraction[i] is area_fractions() at r = 3 * h_tech against the image's
// bounding box when edge_correct, else 1; cells of image -1 keep 1 and get no
// entries. Two phases so the caller can allocate the output exactly:
//   prepare()  -> edge fractions and the column pointer (nnz known),
//   fill()     -> row indices (ascending within each column) and values.
class AmbientFieldBuilder {
 public:
  AmbientFieldBuilder();
  ~AmbientFieldBuilder();

  // x, y, type_code, image: length n, caller-owned, alive until fill() returns.
  // type_code must be >= 0 for every cell with image >= 0.
  Status prepare(Span<const double> x, Span<const double> y, Span<const int> type_code,
                 Span<const int> image, int n_images, double h_tech, bool edge_correct,
                 Span<const double> cos_theta, Span<const double> sin_theta,
                 Span<double> edge_fraction, int n_threads, const InterruptCheck& interrupted);

  std::int64_t nnz() const;
  const std::vector<std::int64_t>& column_pointer() const;

  // row_index and values have length nnz(); edge_fraction is the prepare() output.
  Status fill(Span<int> row_index, Span<double> values, Span<const double> edge_fraction,
              int n_threads, const InterruptCheck& interrupted);

 private:
  struct State;
  std::unique_ptr<State> state_;
};

// Fraction of each cell's neighbours (within `radius`, same image) that share its
// type code, for images with at least `min_image_cells` cells. A cell with no
// neighbours gets 1. Cells in smaller images, or image -1, get NaN.
// The ratio is (long double) same / count, rounded to double, as R's mean() of a
// logical vector.
Status same_type_fraction(Span<const double> x, Span<const double> y,
                          Span<const int> type_code, Span<const int> image,
                          double radius, int min_image_cells, Span<double> fraction,
                          int n_threads, const InterruptCheck& interrupted);

}  // namespace pace

#endif  // PACE_NEIGHBOURHOOD_HPP
