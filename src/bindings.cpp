// bindings.cpp -- thin Rcpp layer over the R-free core in src/core/.
//
// This file only (1) converts R objects to core spans without copying, (2)
// allocates R outputs on the main thread, (3) supplies the interrupt check, and
// (4) turns a core Status into an R error or interrupt. No computation lives here.
#include <Rcpp.h>

#include <cstdint>
#include <string>

#include "core/count_stats.hpp"
#include "core/core_types.hpp"
#include "core/neighbourhood.hpp"

namespace {

// R_CheckUserInterrupt() longjmps; running it under R_ToplevelExec() turns a
// pending interrupt into a FALSE return instead.
void check_interrupt_callback(void*) { R_CheckUserInterrupt(); }

bool user_interrupted() { return R_ToplevelExec(check_interrupt_callback, nullptr) == FALSE; }

// Raise the core status in R: an interrupt stays an interrupt, anything else is an error.
void raise_if_failed(const pace::Status& status, const char* where) {
  if (status.is_ok()) return;
  if (status.code == pace::StatusCode::interrupted) throw Rcpp::internal::InterruptedException();
  Rcpp::stop(std::string(where) + ": " + status.message);
}

pace::Span<const double> column_span(const Rcpp::NumericMatrix& matrix, int column) {
  const std::int64_t n_rows = matrix.nrow();
  return pace::Span<const double>(REAL(matrix) + column * n_rows, n_rows);
}

pace::Span<const int> int_span(const Rcpp::IntegerVector& vector) {
  return pace::Span<const int>(INTEGER(vector), Rf_xlength(vector));
}

pace::Span<const double> double_span(const Rcpp::NumericVector& vector) {
  return pace::Span<const double>(REAL(vector), Rf_xlength(vector));
}

void check_coordinates(const Rcpp::NumericMatrix& coords) {
  if (coords.ncol() != 2) Rcpp::stop("coords must have exactly two columns");
}

}  // namespace

// Neighbour kernels (see pace::neighbour_kernels). coords: n x 2 double;
// neighbour_type, group: 0-based integer codes, -1 for none. Returns a list with
// K_bio and K_tech, n x n_types double matrices.
// [[Rcpp::export]]
Rcpp::List pace_neighbour_kernels_cpp(const Rcpp::NumericMatrix& coords,
                                      const Rcpp::IntegerVector& neighbour_type, int n_types,
                                      const Rcpp::IntegerVector& group, bool per_group,
                                      double h_bio, double h_tech, double eps, int n_threads) {
  check_coordinates(coords);
  const std::int64_t n = coords.nrow();
  Rcpp::NumericMatrix k_bio(n, n_types);
  Rcpp::NumericMatrix k_tech(n, n_types);
  const pace::Status status = pace::neighbour_kernels(
      column_span(coords, 0), column_span(coords, 1), int_span(neighbour_type), n_types,
      int_span(group), per_group, h_bio, h_tech, eps,
      pace::Span<double>(REAL(k_bio), n * n_types), pace::Span<double>(REAL(k_tech), n * n_types),
      n_threads, user_interrupted);
  raise_if_failed(status, "neighbour kernels");
  return Rcpp::List::create(Rcpp::Named("K_bio") = k_bio, Rcpp::Named("K_tech") = k_tech);
}

// Neighbour count per cell under the kernel search rules (test helper).
// [[Rcpp::export]]
Rcpp::IntegerVector pace_neighbour_counts_cpp(const Rcpp::NumericMatrix& coords,
                                              const Rcpp::IntegerVector& group, bool per_group,
                                              double eps, int n_threads) {
  check_coordinates(coords);
  Rcpp::IntegerVector counts(coords.nrow());
  const pace::Status status = pace::neighbour_counts(
      column_span(coords, 0), column_span(coords, 1), int_span(group), per_group, eps,
      pace::Span<int>(INTEGER(counts), Rf_xlength(counts)), n_threads, user_interrupted);
  raise_if_failed(status, "neighbour counts");
  return counts;
}

// Neighbour lists (test helper): list(offsets = n + 1 offsets stored as doubles
// (exact up to 2^53, beyond the int range), neighbours = 1-based cell indices sorted
// within each cell, distances).
// [[Rcpp::export]]
Rcpp::List pace_neighbour_lists_cpp(const Rcpp::NumericMatrix& coords, const Rcpp::IntegerVector& group,
                                    bool per_group, double eps, int n_threads) {
  check_coordinates(coords);
  std::vector<std::int64_t> offsets;
  std::vector<int> neighbours;
  std::vector<double> distances;
  const pace::Status status = pace::neighbour_lists(
      column_span(coords, 0), column_span(coords, 1), int_span(group), per_group, eps,
      offsets, neighbours, distances, n_threads, user_interrupted);
  raise_if_failed(status, "neighbour lists");
  Rcpp::NumericVector offsets_out(offsets.begin(), offsets.end());
  Rcpp::IntegerVector neighbours_out(neighbours.size());
  for (std::size_t k = 0; k < neighbours.size(); ++k) neighbours_out[k] = neighbours[k] + 1;
  Rcpp::NumericVector distances_out(distances.begin(), distances.end());
  return Rcpp::List::create(Rcpp::Named("offsets") = offsets_out,
                            Rcpp::Named("neighbours") = neighbours_out,
                            Rcpp::Named("distances") = distances_out);
}

// Edge fractions for all rows of coords against one rectangle (see pace::area_fractions).
// [[Rcpp::export]]
Rcpp::NumericVector pace_area_fraction_cpp(const Rcpp::NumericMatrix& coords, double r,
                                           double x_min, double x_max, double y_min, double y_max,
                                           const Rcpp::NumericVector& cos_theta,
                                           const Rcpp::NumericVector& sin_theta, int n_threads) {
  check_coordinates(coords);
  const std::int64_t n = coords.nrow();
  Rcpp::NumericVector fraction(n);
  if (n == 0) return fraction;
  Rcpp::IntegerVector rows = Rcpp::seq(0, static_cast<int>(n) - 1);
  const pace::Status status = pace::area_fractions(
      column_span(coords, 0), column_span(coords, 1), int_span(rows), r, x_min, x_max, y_min, y_max,
      double_span(cos_theta), double_span(sin_theta), pace::Span<double>(REAL(fraction), n),
      n_threads, user_interrupted);
  raise_if_failed(status, "area fraction");
  return fraction;
}

// Ambient weight matrix W (see pace::AmbientFieldBuilder). type_code, image: 0-based
// codes (image -1: no image). Returns the dgCMatrix slots i, p, x and the edge fractions.
// [[Rcpp::export]]
Rcpp::List pace_ambient_field_cpp(const Rcpp::NumericMatrix& coords,
                                  const Rcpp::IntegerVector& type_code,
                                  const Rcpp::IntegerVector& image, int n_images, double h_tech,
                                  bool edge_correct, const Rcpp::NumericVector& cos_theta,
                                  const Rcpp::NumericVector& sin_theta, int n_threads) {
  check_coordinates(coords);
  const std::int64_t n = coords.nrow();
  Rcpp::NumericVector edge_fraction(n);
  pace::AmbientFieldBuilder builder;
  pace::Status status = builder.prepare(
      column_span(coords, 0), column_span(coords, 1), int_span(type_code), int_span(image), n_images,
      h_tech, edge_correct, double_span(cos_theta), double_span(sin_theta),
      pace::Span<double>(REAL(edge_fraction), n), n_threads, user_interrupted);
  raise_if_failed(status, "ambient field");

  const std::int64_t nnz = builder.nnz();
  Rcpp::IntegerVector row_index(nnz);
  Rcpp::NumericVector values(nnz);
  Rcpp::IntegerVector column_pointer(n + 1);
  const std::vector<std::int64_t>& pointer = builder.column_pointer();
  for (std::int64_t j = 0; j <= n; ++j) column_pointer[j] = static_cast<int>(pointer[j]);
  status = builder.fill(pace::Span<int>(INTEGER(row_index), nnz), pace::Span<double>(REAL(values), nnz),
                        pace::Span<const double>(REAL(edge_fraction), n), n_threads, user_interrupted);
  raise_if_failed(status, "ambient field");
  return Rcpp::List::create(Rcpp::Named("i") = row_index, Rcpp::Named("p") = column_pointer,
                            Rcpp::Named("x") = values, Rcpp::Named("edge_fraction") = edge_fraction);
}

// Same-type neighbour fraction per cell (see pace::same_type_fraction); NA where
// the cell's image is too small or missing.
// [[Rcpp::export]]
Rcpp::NumericVector pace_same_type_fraction_cpp(const Rcpp::NumericMatrix& coords,
                                                const Rcpp::IntegerVector& type_code,
                                                const Rcpp::IntegerVector& image,
                                                double radius, int min_image_cells, int n_threads) {
  check_coordinates(coords);
  const std::int64_t n = coords.nrow();
  Rcpp::NumericVector fraction(n);
  const pace::Status status = pace::same_type_fraction(
      column_span(coords, 0), column_span(coords, 1), int_span(type_code), int_span(image),
      radius, min_image_cells, pace::Span<double>(REAL(fraction), n), n_threads, user_interrupted);
  raise_if_failed(status, "same-type fraction");
  for (std::int64_t i = 0; i < n; ++i) {
    if (ISNAN(fraction[i])) fraction[i] = NA_REAL;
  }
  return fraction;
}

// Per-group column means of a cells x genes dgCMatrix (see pace::group_column_means).
// Returns an n_groups x n_genes double matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix pace_group_column_means_cpp(const Rcpp::S4& counts, const Rcpp::IntegerVector& group,
                                                int n_groups, bool detection, int n_threads) {
  const Rcpp::IntegerVector dims = counts.slot("Dim");
  const Rcpp::IntegerVector column_pointer = counts.slot("p");
  const Rcpp::IntegerVector row_index = counts.slot("i");
  const Rcpp::NumericVector values = counts.slot("x");
  pace::CscView view;
  view.column_pointer = int_span(column_pointer);
  view.row_index = int_span(row_index);
  view.values = double_span(values);
  view.n_rows = dims[0];
  view.n_cols = dims[1];
  Rcpp::NumericMatrix means(n_groups, dims[1]);
  const pace::Status status = pace::group_column_means(
      view, int_span(group), n_groups, detection,
      pace::Span<double>(REAL(means), static_cast<std::int64_t>(n_groups) * dims[1]), n_threads,
      user_interrupted);
  raise_if_failed(status, "group column means");
  return means;
}
