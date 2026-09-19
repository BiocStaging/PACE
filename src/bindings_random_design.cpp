// The random-effect design builder.
//
// Kept in its own translation unit rather than added to bindings.cpp: it shares
// nothing with the fitting bindings, and Rcpp scans every .cpp in src/ anyway.

#include <Rcpp.h>

#include <cstdint>
#include <vector>

#include "core/random_design.hpp"
#include "core/sparse_product.hpp"
#include "core/count_stats.hpp"

namespace {
// bindings.cpp keeps the same pair in its own anonymous namespace; there is no
// header to share it through. R_CheckUserInterrupt() longjmps, so it runs under
// R_ToplevelExec(), which turns a pending interrupt into a FALSE return.
void check_interrupt_here(void*) { R_CheckUserInterrupt(); }
bool user_interrupted() { return R_ToplevelExec(check_interrupt_here, nullptr) == FALSE; }
}  // namespace

// The sparse random-effect design of one block, as a dgCMatrix, plus the cell
// indices grouped by group.
//
// `x_terms` is the dense model.matrix of the block's formula -- that part stays
// in R, because a formula is an R object and model.matrix() evaluates it -- and
// `cell_group` is a 0-based group code per cell.
//
// What moves here is everything after that: R built n * K_terms triplets in
// three vectors and handed them to Matrix::sparseMatrix(), which sorted them
// into column order. At 1.2M cells and twenty terms those vectors are roughly
// 400 MB and the sort is pure waste, because walking the groups in order and
// the cells ascending within a group emits the rows already sorted. The CSC
// arrays are written once, directly.
//
// [[Rcpp::export]]
Rcpp::List pace_random_design_block_cpp(const Rcpp::NumericMatrix& x_terms,
                                        const Rcpp::IntegerVector& cell_group,
                                        int n_groups) {
  const std::int64_t n = x_terms.nrow();
  const std::int64_t n_terms = x_terms.ncol();
  const std::int64_t n_columns = n_terms * static_cast<std::int64_t>(n_groups);
  const std::int64_t nnz = n * n_terms;

  Rcpp::IntegerVector column_pointer(n_columns + 1);
  Rcpp::IntegerVector row_index(nnz);
  Rcpp::NumericVector values(nnz);
  Rcpp::IntegerVector cells_by_group(n);
  Rcpp::IntegerVector group_start(static_cast<std::int64_t>(n_groups) + 1);

  const pace::Status status = pace::random_design_block(
      pace::Span<const double>(x_terms.begin(), nnz),
      pace::Span<const int>(cell_group.begin(), n),
      n, n_terms, n_groups,
      pace::Span<int>(column_pointer.begin(), n_columns + 1),
      pace::Span<int>(row_index.begin(), nnz),
      pace::Span<double>(values.begin(), nnz),
      pace::Span<int>(cells_by_group.begin(), n),
      pace::Span<int>(group_start.begin(), static_cast<std::int64_t>(n_groups) + 1));
  if (!status.is_ok()) {
    Rcpp::stop("random design block: %s", status.message.c_str());
  }

  Rcpp::S4 z_design("dgCMatrix");
  z_design.slot("i") = row_index;
  z_design.slot("p") = column_pointer;
  z_design.slot("x") = values;
  z_design.slot("Dim") = Rcpp::IntegerVector::create(static_cast<int>(n),
                                                     static_cast<int>(n_columns));
  z_design.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
  z_design.slot("factors") = Rcpp::List::create();

  return Rcpp::List::create(Rcpp::Named("Z") = z_design,
                            Rcpp::Named("cells_by_group") = cells_by_group,
                            Rcpp::Named("group_start") = group_start);
}

