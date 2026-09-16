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
#include "core/decomposition.hpp"
#include "core/hyperparameters.hpp"
#include "core/irls_chunk.hpp"
#include "core/neighbourhood.hpp"
#include "core/statistics.hpp"

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

// The dgCMatrix slots, held so the core's view stays valid for the call.
struct CscHolder {
  Rcpp::IntegerVector dims;
  Rcpp::IntegerVector column_pointer;
  Rcpp::IntegerVector row_index;
  Rcpp::NumericVector values;
  pace::CscView view;

  explicit CscHolder(const Rcpp::S4& matrix)
      : dims(matrix.slot("Dim")),
        column_pointer(matrix.slot("p")),
        row_index(matrix.slot("i")),
        values(matrix.slot("x")) {
    view.column_pointer = int_span(column_pointer);
    view.row_index = int_span(row_index);
    view.values = double_span(values);
    view.n_rows = dims[0];
    view.n_cols = dims[1];
  }
};

pace::Span<double> out_span(Rcpp::NumericVector& vector) {
  return pace::Span<double>(REAL(vector), Rf_xlength(vector));
}

pace::Span<int> out_span(Rcpp::IntegerVector& vector) {
  return pace::Span<int>(INTEGER(vector), Rf_xlength(vector));
}

pace::Span<double> out_span(Rcpp::NumericMatrix& matrix) {
  return pace::Span<double>(REAL(matrix), Rf_xlength(matrix));
}

pace::Span<const double> const_span(const Rcpp::NumericMatrix& matrix) {
  return pace::Span<const double>(REAL(matrix), Rf_xlength(matrix));
}

pace::Span<const int> const_span(const Rcpp::IntegerMatrix& matrix) {
  return pace::Span<const int>(INTEGER(matrix), Rf_xlength(matrix));
}

// R's NA for the slots where stats::var() has fewer than two values.
void set_na_where_too_few(Rcpp::NumericVector& variance, const Rcpp::IntegerVector& n_used) {
  for (R_xlen_t k = 0; k < variance.size(); ++k) {
    if (n_used[k] < 2) variance[k] = NA_REAL;
  }
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

// Per-group covariance of the columns of a dense n x p matrix (see
// pace::group_covariances). Returns a p x p x n_groups array (slice [, , g] is
// group g's covariance matrix), matching the core's contiguous per-group layout.
// [[Rcpp::export]]
Rcpp::NumericVector pace_group_covariances_cpp(const Rcpp::NumericMatrix& values,
                                               const Rcpp::IntegerVector& group, int n_groups,
                                               int n_threads) {
  const std::int64_t n = values.nrow();
  const std::int64_t p = values.ncol();
  Rcpp::NumericVector covariance(static_cast<R_xlen_t>(n_groups) * p * p);
  const pace::Status status = pace::group_covariances(
      pace::Span<const double>(REAL(values), n * p), n, p, int_span(group), n_groups,
      pace::Span<double>(REAL(covariance), Rf_xlength(covariance)), n_threads, user_interrupted);
  raise_if_failed(status, "group covariances");
  covariance.attr("dim") = Rcpp::IntegerVector::create(static_cast<int>(p), static_cast<int>(p), n_groups);
  return covariance;
}

// Single-frame statistics of a cells x genes dgCMatrix (see
// pace::single_frame_statistics). Returns focal_mean and within_ss
// (n_groups x n_genes) and global_mean (n_genes).
// [[Rcpp::export]]
Rcpp::List pace_single_frame_statistics_cpp(const Rcpp::S4& counts, const Rcpp::NumericVector& n_count,
                                            const Rcpp::IntegerVector& group, int n_groups,
                                            int n_threads) {
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
  Rcpp::NumericMatrix focal_mean(n_groups, dims[1]);
  Rcpp::NumericMatrix within_ss(n_groups, dims[1]);
  Rcpp::NumericVector global_mean(dims[1]);
  const std::int64_t n_slots = static_cast<std::int64_t>(n_groups) * dims[1];
  const pace::Status status = pace::single_frame_statistics(
      view, double_span(n_count), int_span(group), n_groups,
      pace::Span<double>(REAL(focal_mean), n_slots), pace::Span<double>(REAL(within_ss), n_slots),
      pace::Span<double>(REAL(global_mean), dims[1]), n_threads, user_interrupted);
  raise_if_failed(status, "single-frame statistics");
  return Rcpp::List::create(Rcpp::Named("focal_mean") = focal_mean,
                            Rcpp::Named("within_ss") = within_ss,
                            Rcpp::Named("global_mean") = global_mean);
}

namespace {

// A gene range of a cells x genes CSC matrix; `first_gene` is 1-based.
pace::GeneBlock gene_block(const CscHolder& holder, int first_gene) {
  pace::GeneBlock block;
  block.matrix = holder.view;
  block.first_gene = first_gene - 1;
  return block;
}

}  // namespace

// Per-group mean and variance of the columns of a dense matrix (see
// pace::dense_group_moments). `values` may be a 0 x 0 matrix, which means an
// all-zero n x p matrix. Returns mean, variance (R's NA where fewer than two
// values entered it) and any_nonzero.
// [[Rcpp::export]]
Rcpp::List pace_dense_group_moments_cpp(const Rcpp::NumericMatrix& values, int n, int p,
                                        const Rcpp::IntegerVector& group, int n_groups,
                                        bool want_mean, bool want_variance, int n_threads) {
  const bool all_zero = values.size() == 0;
  if (!all_zero && (values.nrow() != n || values.ncol() != p)) {
    Rcpp::stop("values must be empty or n x p");
  }
  Rcpp::NumericMatrix mean(want_mean ? n_groups : 0, want_mean ? p : 0);
  Rcpp::NumericMatrix variance(want_variance ? n_groups : 0, want_variance ? p : 0);
  Rcpp::IntegerVector n_used(want_variance ? static_cast<R_xlen_t>(n_groups) * p : 0);
  bool any_nonzero = false;
  const pace::Status status = pace::dense_group_moments(
      all_zero ? pace::Span<const double>() : const_span(values), n, p, int_span(group), n_groups,
      out_span(mean), out_span(variance), out_span(n_used), &any_nonzero, n_threads, user_interrupted);
  raise_if_failed(status, "group moments");
  if (want_variance) set_na_where_too_few(variance, n_used);
  return Rcpp::List::create(Rcpp::Named("mean") = mean, Rcpp::Named("variance") = variance,
                            Rcpp::Named("any_nonzero") = any_nonzero);
}

// One gene chunk of the streaming solver's final pass (see pace::final_pass_statistics).
// `return_matrices` also returns this chunk's mu and technical offset.
// [[Rcpp::export]]
Rcpp::List pace_final_pass_statistics_cpp(const Rcpp::NumericMatrix& eta,
                                          const Rcpp::NumericVector& offset,
                                          const Rcpp::S4& ambient, int first_gene,
                                          const Rcpp::NumericVector& rho,
                                          const Rcpp::IntegerVector& group, int n_groups,
                                          bool want_groups, bool return_matrices) {
  CscHolder ambient_holder(ambient);
  const std::int64_t n = eta.nrow();
  const std::int64_t n_genes = eta.ncol();
  Rcpp::NumericMatrix mu_group_sum(want_groups ? n_groups : 0, want_groups ? n_genes : 0);
  Rcpp::NumericMatrix toff_variance(want_groups ? n_groups : 0, want_groups ? n_genes : 0);
  Rcpp::NumericVector mu_column_sum(n_genes);
  Rcpp::NumericVector spill_row_sum(n);
  Rcpp::NumericVector total_row_sum(n);
  Rcpp::NumericMatrix mu(return_matrices ? n : 0, return_matrices ? n_genes : 0);
  Rcpp::NumericMatrix toff(return_matrices ? n : 0, return_matrices ? n_genes : 0);
  bool any_nonzero = false;
  const pace::Status status = pace::final_pass_statistics(
      const_span(eta), n, n_genes, double_span(offset), gene_block(ambient_holder, first_gene),
      double_span(rho), int_span(group), n_groups, out_span(mu_group_sum), out_span(toff_variance),
      out_span(mu_column_sum), out_span(spill_row_sum), out_span(total_row_sum), out_span(mu),
      out_span(toff), &any_nonzero, user_interrupted);
  raise_if_failed(status, "final pass statistics");
  // stats::var() of fewer than two values is NA, not NaN.
  std::vector<int> group_size(want_groups ? n_groups : 0, 0);
  if (want_groups) {
    for (R_xlen_t i = 0; i < group.size(); ++i) {
      const int g = group[i];
      if (g >= 0 && g < n_groups) group_size[g] += 1;
    }
    for (std::int64_t j = 0; j < n_genes; ++j) {
      for (int g = 0; g < n_groups; ++g) {
        if (group_size[g] < 2) toff_variance[g + j * n_groups] = NA_REAL;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("mu_group_sum") = mu_group_sum,
                            Rcpp::Named("toff_variance") = toff_variance,
                            Rcpp::Named("mu_column_sum") = mu_column_sum,
                            Rcpp::Named("spill_row_sum") = spill_row_sum,
                            Rcpp::Named("total_row_sum") = total_row_sum,
                            Rcpp::Named("any_nonzero") = any_nonzero,
                            Rcpp::Named("mu") = mu, Rcpp::Named("toff") = toff);
}

// Per-(focal, gene) variance decomposition (see pace::variance_decomposition).
// Index arguments are 1-based; a 0 x 0 matrix means the block is absent. Every
// returned vector has n_focals * n_genes entries, focal-major.
// [[Rcpp::export]]
Rcpp::List pace_variance_decomposition_cpp(
    const Rcpp::NumericMatrix& ct_means, const Rcpp::IntegerVector& group_size,
    const Rcpp::IntegerVector& focal_group, const Rcpp::IntegerVector& n_focal,
    const Rcpp::NumericMatrix& u, const Rcpp::NumericMatrix& se_u,
    const Rcpp::IntegerMatrix& slope_rows, const Rcpp::NumericVector& kernel_cov,
    const Rcpp::IntegerMatrix& responder_rows, const Rcpp::IntegerVector& responder_keep,
    const Rcpp::NumericVector& responder_cov, const Rcpp::IntegerVector& intercept_rows,
    const Rcpp::NumericMatrix& toff_var, const Rcpp::NumericVector& spill_cov,
    const Rcpp::NumericMatrix& beta_spill, const Rcpp::NumericMatrix& mu_mean,
    const Rcpp::NumericVector& alpha, bool nb1, int n_threads) {
  pace::DecompositionInput input;
  input.n_genes = ct_means.ncol();
  input.n_groups = group_size.size();
  input.n_focals = focal_group.size();
  input.ct_means = const_span(ct_means);
  input.group_size = int_span(group_size);
  input.focal_group = int_span(focal_group);
  input.n_focal = int_span(n_focal);
  input.u = const_span(u);
  input.se_u = const_span(se_u);
  input.q = u.nrow();
  input.n_kernel = slope_rows.nrow();
  input.slope_rows = const_span(slope_rows);
  input.kernel_cov = double_span(kernel_cov);
  input.n_responder = responder_rows.nrow();
  input.responder_rows = const_span(responder_rows);
  input.responder_keep = int_span(responder_keep);
  input.responder_cov = double_span(responder_cov);
  input.intercept_rows = int_span(intercept_rows);
  input.toff_var = const_span(toff_var);
  input.n_spill = beta_spill.nrow();
  input.spill_cov = double_span(spill_cov);
  input.beta_spill = const_span(beta_spill);
  input.mu_mean = const_span(mu_mean);
  input.alpha = double_span(alpha);
  input.nb1 = nb1;

  const R_xlen_t rows = static_cast<R_xlen_t>(input.n_focals) * input.n_genes;
  Rcpp::NumericVector v_state_baseline(rows), v_state_responder(rows), v_spill(rows), v_disp(rows);
  Rcpp::NumericVector celltype_offset_sq(rows), focal_mean(rows), max_other_mean(rows), spec(rows);
  Rcpp::NumericVector focal_other_ratio(rows), pct_celltype(rows), pct_state(rows);
  Rcpp::NumericVector pct_responder(rows), pct_spill(rows), pct_residual(rows);
  Rcpp::IntegerVector is_contaminated(rows), keep(rows);
  pace::DecompositionOutput output;
  output.v_state_baseline = out_span(v_state_baseline);
  output.v_state_responder = out_span(v_state_responder);
  output.v_spill = out_span(v_spill);
  output.v_disp = out_span(v_disp);
  output.celltype_offset_sq = out_span(celltype_offset_sq);
  output.focal_mean = out_span(focal_mean);
  output.max_other_mean = out_span(max_other_mean);
  output.spec = out_span(spec);
  output.focal_other_ratio = out_span(focal_other_ratio);
  output.is_contaminated = out_span(is_contaminated);
  output.pct_celltype = out_span(pct_celltype);
  output.pct_state = out_span(pct_state);
  output.pct_responder = out_span(pct_responder);
  output.pct_spill = out_span(pct_spill);
  output.pct_residual = out_span(pct_residual);
  output.keep = out_span(keep);
  const pace::Status status = pace::variance_decomposition(input, output, n_threads, user_interrupted);
  raise_if_failed(status, "variance decomposition");
  return Rcpp::List::create(
      Rcpp::Named("V_state_baseline") = v_state_baseline,
      Rcpp::Named("V_state_responder") = v_state_responder, Rcpp::Named("V_spill") = v_spill,
      Rcpp::Named("V_disp") = v_disp, Rcpp::Named("celltype_offset_sq") = celltype_offset_sq,
      Rcpp::Named("focal_mean") = focal_mean, Rcpp::Named("max_other_mean") = max_other_mean,
      Rcpp::Named("spec") = spec, Rcpp::Named("focal_other_ratio") = focal_other_ratio,
      Rcpp::Named("is_contaminated") = is_contaminated, Rcpp::Named("pct_celltype") = pct_celltype,
      Rcpp::Named("pct_state") = pct_state, Rcpp::Named("pct_responder") = pct_responder,
      Rcpp::Named("pct_spill") = pct_spill, Rcpp::Named("pct_residual") = pct_residual,
      Rcpp::Named("keep") = keep);
}

// The aggregate tables of the decomposition (see pace::decomposition_aggregates).
// `focal_code` is 0-based, in the row order the caller wants the tables in.
// [[Rcpp::export]]
Rcpp::List pace_decomposition_aggregates_cpp(const Rcpp::IntegerVector& focal_code, int n_focals,
                                             const Rcpp::NumericMatrix& pct,
                                             const Rcpp::NumericMatrix& components,
                                             const Rcpp::NumericVector& spec,
                                             const Rcpp::IntegerVector& is_contaminated,
                                             bool weight_by_spec_sq) {
  const std::int64_t n_rows = focal_code.size();
  Rcpp::NumericMatrix mean_table(n_focals, 5);
  Rcpp::NumericMatrix specw_table(weight_by_spec_sq ? n_focals : 0, weight_by_spec_sq ? 5 : 0);
  Rcpp::NumericMatrix pooled_table(n_focals, 5);
  Rcpp::NumericVector total(n_rows);
  Rcpp::NumericVector total_ss(n_focals);
  Rcpp::IntegerVector n_genes(n_focals), n_specific(n_focals);
  const pace::Status status = pace::decomposition_aggregates(
      int_span(focal_code), n_focals, n_rows, const_span(pct), const_span(components),
      double_span(spec), int_span(is_contaminated), out_span(mean_table), out_span(specw_table),
      out_span(pooled_table), out_span(total), out_span(total_ss), out_span(n_genes),
      out_span(n_specific));
  raise_if_failed(status, "decomposition aggregates");
  return Rcpp::List::create(Rcpp::Named("mean") = mean_table, Rcpp::Named("specw") = specw_table,
                            Rcpp::Named("pooled") = pooled_table, Rcpp::Named("Total") = total,
                            Rcpp::Named("total_SS") = total_ss, Rcpp::Named("n_genes") = n_genes,
                            Rcpp::Named("n_specific_genes") = n_specific);
}

// The four-block view of the decomposition (see pace::four_block_shares).
// [[Rcpp::export]]
Rcpp::List pace_four_block_shares_cpp(const Rcpp::NumericVector& celltype_offset_sq,
                                      const Rcpp::NumericVector& v_state_baseline,
                                      const Rcpp::NumericVector& v_state_responder,
                                      const Rcpp::NumericVector& v_spill,
                                      const Rcpp::NumericVector& v_disp) {
  const std::int64_t n_rows = celltype_offset_sq.size();
  Rcpp::NumericVector total(n_rows), v_state(n_rows), pct_celltype(n_rows), pct_state(n_rows);
  Rcpp::NumericVector pct_spill(n_rows), pct_residual(n_rows);
  const pace::Status status = pace::four_block_shares(
      double_span(celltype_offset_sq), double_span(v_state_baseline), double_span(v_state_responder),
      double_span(v_spill), double_span(v_disp), n_rows, out_span(total), out_span(v_state),
      out_span(pct_celltype), out_span(pct_state), out_span(pct_spill), out_span(pct_residual));
  raise_if_failed(status, "four-block shares");
  return Rcpp::List::create(Rcpp::Named("Total") = total, Rcpp::Named("V_state") = v_state,
                            Rcpp::Named("pct_celltype") = pct_celltype,
                            Rcpp::Named("pct_state") = pct_state,
                            Rcpp::Named("pct_spill") = pct_spill,
                            Rcpp::Named("pct_residual") = pct_residual);
}

// The observed single-frame shares (see pace::single_frame_shares). The four V
// matrices are n_groups x n_genes, aligned to the same genes as the statistics.
// Returned vectors are n_groups * n_genes, group-major.
// [[Rcpp::export]]
Rcpp::List pace_single_frame_shares_cpp(const Rcpp::NumericMatrix& focal_mean,
                                        const Rcpp::NumericMatrix& within_ss,
                                        const Rcpp::NumericVector& global_mean,
                                        const Rcpp::IntegerVector& group_size,
                                        const Rcpp::NumericMatrix& v_state,
                                        const Rcpp::NumericMatrix& v_responder,
                                        const Rcpp::NumericMatrix& v_spill,
                                        const Rcpp::NumericMatrix& v_disp, int n_threads) {
  const int n_groups = group_size.size();
  const std::int64_t n_genes = focal_mean.ncol();
  const R_xlen_t rows = static_cast<R_xlen_t>(n_groups) * n_genes;
  const bool has_responder = v_responder.size() != 0;
  Rcpp::NumericVector pct_celltype(rows), pct_spatial(rows), pct_spill(rows), pct_residual(rows);
  Rcpp::NumericVector pct_responder(has_responder ? rows : 0);
  Rcpp::NumericVector ss_lineage(rows), ss_within(rows), denom(rows);
  Rcpp::IntegerVector keep(n_groups);
  const pace::Status status = pace::single_frame_shares(
      const_span(focal_mean), const_span(within_ss), double_span(global_mean), int_span(group_size),
      n_groups, n_genes, const_span(v_state), const_span(v_responder), const_span(v_spill),
      const_span(v_disp), out_span(pct_celltype), out_span(pct_spatial), out_span(pct_responder),
      out_span(pct_spill), out_span(pct_residual), out_span(ss_lineage), out_span(ss_within),
      out_span(denom), out_span(keep), n_threads, user_interrupted);
  raise_if_failed(status, "single-frame shares");
  return Rcpp::List::create(Rcpp::Named("pct_celltype") = pct_celltype,
                            Rcpp::Named("pct_spatial") = pct_spatial,
                            Rcpp::Named("pct_responder") = pct_responder,
                            Rcpp::Named("pct_spill") = pct_spill,
                            Rcpp::Named("pct_residual") = pct_residual,
                            Rcpp::Named("SS_lineage") = ss_lineage,
                            Rcpp::Named("SS_within") = ss_within, Rcpp::Named("denom") = denom,
                            Rcpp::Named("keep") = keep);
}

// Driver scores for one (focal, neighbour) pair (see pace::driver_scores).
// [[Rcpp::export]]
Rcpp::List pace_driver_scores_cpp(const Rcpp::NumericVector& estimate_shrunk,
                                  const Rcpp::NumericVector& spec,
                                  const Rcpp::NumericVector& focal_mean,
                                  const Rcpp::NumericVector& mu_bar,
                                  const Rcpp::NumericVector& alpha,
                                  const Rcpp::NumericVector& u_raw, double var_n, double var_rn,
                                  bool has_responder) {
  const std::int64_t n = estimate_shrunk.size();
  Rcpp::NumericVector mcsd(n), mcsd4(n), v_resid(n), v_s(n), v_total(n), r2_s(n);
  Rcpp::NumericVector v_rxs(has_responder ? n : 0), r2_rxs(has_responder ? n : 0);
  const pace::Status status = pace::driver_scores(
      double_span(estimate_shrunk), double_span(spec), double_span(focal_mean), double_span(mu_bar),
      double_span(alpha), double_span(u_raw), var_n, var_rn, has_responder, n, out_span(mcsd),
      out_span(mcsd4), out_span(v_resid), out_span(v_s), out_span(v_rxs), out_span(v_total),
      out_span(r2_s), out_span(r2_rxs));
  raise_if_failed(status, "driver scores");
  return Rcpp::List::create(Rcpp::Named("MCSD") = mcsd, Rcpp::Named("MCSD4") = mcsd4,
                            Rcpp::Named("V_resid") = v_resid, Rcpp::Named("V_S") = v_s,
                            Rcpp::Named("V_RxS") = v_rxs, Rcpp::Named("V_total") = v_total,
                            Rcpp::Named("R2_S") = r2_s, Rcpp::Named("R2_RxS") = r2_rxs);
}

// stats::cov() of selected columns of a sparse matrix over selected rows
// (see pace::subset_covariance). `rows` and `cols` are 1-based.
// [[Rcpp::export]]
Rcpp::NumericMatrix pace_subset_covariance_cpp(const Rcpp::S4& matrix, const Rcpp::IntegerVector& rows,
                                               const Rcpp::IntegerVector& cols,
                                               const Rcpp::NumericVector& scale) {
  CscHolder holder(matrix);
  Rcpp::NumericMatrix covariance(cols.size(), cols.size());
  const pace::Status status = pace::subset_covariance(holder.view, int_span(rows), int_span(cols),
                                                      double_span(scale), out_span(covariance));
  raise_if_failed(status, "subset covariance");
  return covariance;
}

// The values of one (1-based) column of a sparse matrix at selected (1-based) rows.
// [[Rcpp::export]]
Rcpp::NumericVector pace_subset_column_cpp(const Rcpp::S4& matrix, const Rcpp::IntegerVector& rows,
                                           int column) {
  CscHolder holder(matrix);
  Rcpp::NumericVector values(rows.size());
  const pace::Status status = pace::subset_column(holder.view, int_span(rows), column, out_span(values));
  raise_if_failed(status, "subset column");
  return values;
}

// The 1-based rows at which one column of a sparse matrix is not zero.
// [[Rcpp::export]]
Rcpp::IntegerVector pace_column_nonzero_rows_cpp(const Rcpp::S4& matrix, int column) {
  CscHolder holder(matrix);
  std::int64_t n_found = 0;
  pace::Status status = pace::column_nonzero_rows(holder.view, column, pace::Span<int>(), &n_found);
  raise_if_failed(status, "column non-zero rows");
  Rcpp::IntegerVector rows(n_found);
  status = pace::column_nonzero_rows(holder.view, column, out_span(rows), &n_found);
  raise_if_failed(status, "column non-zero rows");
  return rows;
}

// Pratt's pair attribution for one focal cell type (see pace::pair_variance_pratt).
// [[Rcpp::export]]
Rcpp::List pace_pair_variance_pratt_cpp(const Rcpp::NumericMatrix& sigma,
                                        const Rcpp::NumericMatrix& u) {
  const std::int64_t n_pairs = sigma.nrow();
  const std::int64_t n_genes = u.ncol();
  Rcpp::NumericVector v_pair(n_pairs), v_diag(n_pairs), share_pct(n_pairs), diag_pct(n_pairs);
  Rcpp::NumericVector sign(n_pairs);
  double v_total = 0;
  double v_diag_total = 0;
  double cross_cov_pct = 0;
  int n_negative = 0;
  const pace::Status status = pace::pair_variance_pratt(
      const_span(sigma), const_span(u), n_pairs, n_genes, out_span(v_pair), out_span(v_diag),
      out_span(share_pct), out_span(diag_pct), out_span(sign), &v_total, &v_diag_total,
      &cross_cov_pct, &n_negative);
  raise_if_failed(status, "pair variance");
  return Rcpp::List::create(Rcpp::Named("V_pair_pratt") = v_pair, Rcpp::Named("V_pair_diag") = v_diag,
                            Rcpp::Named("within_focal_share_pct") = share_pct,
                            Rcpp::Named("within_focal_diag_pct") = diag_pct,
                            Rcpp::Named("sign") = sign, Rcpp::Named("V_block_pratt") = v_total,
                            Rcpp::Named("V_block_diag") = v_diag_total,
                            Rcpp::Named("cross_cov_pct") = cross_cov_pct,
                            Rcpp::Named("n_negative_pairs") =
                                n_negative < 0 ? NA_INTEGER : n_negative);
}

// SS-weighted per-focal blocks of the single frame (see pace::single_frame_focal_blocks).
// [[Rcpp::export]]
Rcpp::NumericMatrix pace_single_frame_focal_blocks_cpp(const Rcpp::IntegerVector& focal_code,
                                                       int n_focals, const Rcpp::NumericMatrix& pct,
                                                       const Rcpp::NumericVector& denom) {
  Rcpp::NumericMatrix table(n_focals, pct.ncol());
  const pace::Status status = pace::single_frame_focal_blocks(
      int_span(focal_code), n_focals, focal_code.size(), const_span(pct), pct.ncol(),
      double_span(denom), out_span(table));
  raise_if_failed(status, "single-frame focal blocks");
  return table;
}

// True when every value is finite (see pace::all_finite): the same check as
// all(is.finite(x)) without allocating the logical vector.
// [[Rcpp::export]]
bool pace_all_finite_cpp(const Rcpp::NumericVector& values) {
  return pace::all_finite(double_span(values));
}

// ---------------------------------------------------------------------------
// The streaming solver's per-chunk arithmetic (see core/irls_chunk.hpp) and its
// variance-component updates (see core/hyperparameters.hpp). `first_gene` is the
// 1-based first column of the chunk in the counts / ambient matrices.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List pace_working_response_cpp(const Rcpp::NumericMatrix& eta, const Rcpp::S4& counts,
                                     const Rcpp::S4& ambient, int first_gene, int n_genes,
                                     const Rcpp::NumericVector& offset,
                                     const Rcpp::NumericVector& rho,
                                     const Rcpp::NumericVector& alpha,
                                     const Rcpp::NumericVector& sample_weight, bool nb2,
                                     bool seed_iteration, int n_cells, int n_threads) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  Rcpp::NumericMatrix z(n_cells, n_genes);
  Rcpp::NumericMatrix w(n_cells, n_genes);
  Rcpp::NumericVector colsum_w(n_genes);
  const pace::Status status = pace::working_response(
      seed_iteration ? pace::Span<const double>() : const_span(eta), gene_block(count_holder, first_gene),
      gene_block(ambient_holder, first_gene), double_span(offset), double_span(rho),
      double_span(alpha), double_span(sample_weight), nb2, seed_iteration, n_cells, n_genes,
      out_span(z), out_span(w), out_span(colsum_w), n_threads, user_interrupted);
  raise_if_failed(status, "working response");
  return Rcpp::List::create(Rcpp::Named("z") = z, Rcpp::Named("w") = w,
                            Rcpp::Named("colsum_w") = colsum_w);
}

// [[Rcpp::export]]
Rcpp::List pace_rho_accumulate_cpp(const Rcpp::NumericMatrix& eta,
                                   const Rcpp::NumericMatrix& prev_eta, const Rcpp::S4& counts,
                                   const Rcpp::S4& ambient, int first_gene,
                                   const Rcpp::NumericVector& offset,
                                   const Rcpp::NumericVector& rho,
                                   const Rcpp::NumericVector& alpha,
                                   const Rcpp::NumericMatrix& mask,
                                   const Rcpp::IntegerVector& mask_index,
                                   const Rcpp::NumericVector& num_in,
                                   const Rcpp::NumericVector& den_in, bool nb2,
                                   bool seed_iteration, bool seed_previous, int n_cells,
                                   int n_genes_in, bool want_tail_counts) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  const std::int64_t n = seed_iteration ? n_cells : eta.nrow();
  const std::int64_t n_genes = seed_iteration ? n_genes_in : eta.ncol();
  Rcpp::NumericVector num(Rcpp::clone(num_in));
  Rcpp::NumericVector den(Rcpp::clone(den_in));
  double rel_delta_max = 0;
  double rel_delta_sum = 0;
  std::int64_t n_finite = 0;
  std::int64_t n_nonfinite = 0;
  Rcpp::NumericVector tail_counts(want_tail_counts ? 4 : 0);
  const pace::Status status = pace::rho_accumulate(
      const_span(eta), const_span(prev_eta), gene_block(count_holder, first_gene),
      gene_block(ambient_holder, first_gene), double_span(offset), double_span(rho),
      double_span(alpha), const_span(mask), int_span(mask_index), mask.nrow(), nb2, seed_iteration,
      seed_previous, n, n_genes, out_span(num), out_span(den), &rel_delta_max, &rel_delta_sum, &n_finite,
      &n_nonfinite, out_span(tail_counts), user_interrupted);
  raise_if_failed(status, "rho accumulation");
  return Rcpp::List::create(Rcpp::Named("num") = num, Rcpp::Named("den") = den,
                            Rcpp::Named("rel_delta_max") = rel_delta_max,
                            Rcpp::Named("rel_delta_sum") = rel_delta_sum,
                            Rcpp::Named("n_finite") = static_cast<double>(n_finite),
                            Rcpp::Named("n_nonfinite") = static_cast<double>(n_nonfinite),
                            Rcpp::Named("tail_counts") = tail_counts);
}

// [[Rcpp::export]]
Rcpp::List pace_fitted_mean_column_cpp(const Rcpp::NumericVector& eta_column,
                                       const Rcpp::S4& counts, const Rcpp::S4& ambient,
                                       int first_gene, int gene,
                                       const Rcpp::NumericVector& offset,
                                       const Rcpp::NumericVector& rho) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  const std::int64_t n = eta_column.size();
  Rcpp::NumericVector mu(n);
  Rcpp::NumericVector y(n);
  const pace::Status status = pace::fitted_mean_column(
      double_span(eta_column), gene_block(count_holder, first_gene),
      gene_block(ambient_holder, first_gene), double_span(offset), double_span(rho), n, gene - 1,
      out_span(mu), out_span(y));
  raise_if_failed(status, "fitted mean column");
  return Rcpp::List::create(Rcpp::Named("mu") = mu, Rcpp::Named("y") = y);
}

// [[Rcpp::export]]
Rcpp::List pace_rho_shrink_cpp(const Rcpp::NumericVector& num, const Rcpp::NumericVector& den) {
  const std::int64_t n = num.size();
  Rcpp::NumericVector rho(n);
  std::int64_t n_nonfinite = 0;
  double rho_bar = 0;
  double den_quantile = 0;
  const pace::Status status = pace::rho_shrink(double_span(num), double_span(den), n,
                                               out_span(rho), &n_nonfinite, &rho_bar, &den_quantile);
  raise_if_failed(status, "rho shrink");
  return Rcpp::List::create(Rcpp::Named("rho") = rho,
                            Rcpp::Named("n_nonfinite") = static_cast<double>(n_nonfinite),
                            Rcpp::Named("rho_bar") = rho_bar,
                            Rcpp::Named("den_quantile") = den_quantile);
}

// [[Rcpp::export]]
Rcpp::NumericVector pace_tau_em_update_cpp(const Rcpp::NumericMatrix& u,
                                           const Rcpp::NumericMatrix& re_var) {
  Rcpp::NumericVector tau(u.nrow());
  const pace::Status status = pace::tau_em_update(const_span(u), const_span(re_var), u.nrow(),
                                                  u.ncol(), out_span(tau));
  raise_if_failed(status, "tau EM update");
  return tau;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pace_tau_hierarchical_cpp(const Rcpp::NumericMatrix& tau,
                                              const Rcpp::IntegerVector& n_group,
                                              double lambda_factor) {
  Rcpp::NumericMatrix out(tau.nrow(), tau.ncol());
  const pace::Status status = pace::tau_hierarchical(const_span(tau), int_span(n_group), tau.nrow(),
                                                     tau.ncol(), lambda_factor, out_span(out));
  raise_if_failed(status, "hierarchical tau");
  return out;
}

// [[Rcpp::export]]
Rcpp::List pace_tau_eb_summaries_cpp(const Rcpp::NumericMatrix& s2, double reml_factor) {
  const std::int64_t q = s2.nrow();
  Rcpp::NumericVector panel(q), panel_median(q), log_variance(q);
  Rcpp::IntegerVector n_log_finite(q);
  const pace::Status status = pace::tau_eb_summaries(const_span(s2), q, s2.ncol(), reml_factor,
                                                     out_span(panel), out_span(panel_median),
                                                     out_span(log_variance), out_span(n_log_finite));
  raise_if_failed(status, "tau EB summaries");
  return Rcpp::List::create(Rcpp::Named("panel") = panel,
                            Rcpp::Named("panel_median") = panel_median,
                            Rcpp::Named("log_variance") = log_variance,
                            Rcpp::Named("n_log_finite") = n_log_finite);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pace_tau_eb_apply_cpp(const Rcpp::NumericMatrix& s2, double reml_factor,
                                          const Rcpp::NumericVector& panel,
                                          const Rcpp::NumericVector& panel_median,
                                          const Rcpp::NumericVector& d0, double tau_floor) {
  Rcpp::NumericMatrix out(s2.nrow(), s2.ncol());
  const pace::Status status = pace::tau_eb_apply(const_span(s2), s2.nrow(), s2.ncol(), reml_factor,
                                                 double_span(panel), double_span(panel_median),
                                                 double_span(d0), tau_floor, out_span(out));
  raise_if_failed(status, "tau EB shrinkage");
  return out;
}

// [[Rcpp::export]]
Rcpp::List pace_tau_half_cauchy_cpp(const Rcpp::NumericMatrix& s2, int n_em_iter, double tau_floor,
                                    const Rcpp::NumericVector& lambda2_prev,
                                    const Rcpp::NumericMatrix& a_prev) {
  Rcpp::NumericMatrix out(s2.nrow(), s2.ncol());
  Rcpp::NumericMatrix a_out(s2.nrow(), s2.ncol());
  Rcpp::NumericVector lambda2_out(s2.nrow());
  Rcpp::NumericVector panel(s2.nrow());
  const pace::Status status = pace::tau_half_cauchy(
      const_span(s2), s2.nrow(), s2.ncol(), n_em_iter, tau_floor, double_span(lambda2_prev),
      const_span(a_prev), out_span(out), out_span(lambda2_out), out_span(a_out), out_span(panel));
  raise_if_failed(status, "half-Cauchy tau");
  return Rcpp::List::create(Rcpp::Named("tau") = out, Rcpp::Named("lambda_sq") = lambda2_out,
                            Rcpp::Named("a") = a_out, Rcpp::Named("panel") = panel);
}

// [[Rcpp::export]]
Rcpp::List pace_tau_clamp_cpp(const Rcpp::NumericMatrix& tau, double tau_max) {
  Rcpp::NumericMatrix out(Rcpp::clone(tau));
  std::int64_t n_binding = 0;
  const pace::Status status = pace::tau_clamp(out_span(out), Rf_xlength(out), tau_max, &n_binding);
  raise_if_failed(status, "tau clamp");
  return Rcpp::List::create(Rcpp::Named("tau") = out,
                            Rcpp::Named("n_binding") = static_cast<double>(n_binding));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pace_data_informed_weights_cpp(const Rcpp::NumericMatrix& detection_rate,
                                                   const Rcpp::IntegerVector& focal_of_row,
                                                   const Rcpp::NumericVector& scale_of_row) {
  const std::int64_t q = focal_of_row.size();
  const std::int64_t n_genes = detection_rate.ncol();
  Rcpp::NumericMatrix weights(q, n_genes);
  const pace::Status status = pace::data_informed_weights(
      const_span(detection_rate), detection_rate.nrow(), n_genes, int_span(focal_of_row),
      double_span(scale_of_row), q, out_span(weights));
  raise_if_failed(status, "data-informed weights");
  return weights;
}

