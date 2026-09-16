// preprocess.hpp -- the arithmetic between the neighbour kernels and the solver:
// dropping kernel columns that carry no usable variation, centring them within
// each (image, cell type), standardising them for the image block, and the
// anchor decision. These were R loops over cell types; the formulas are
// unchanged, including which mean each one uses.
#ifndef PACE_PREPROCESS_HPP
#define PACE_PREPROCESS_HPP

#include <cstdint>

#include "core_types.hpp"

namespace pace {

// Zero the (focal, neighbour) kernel columns with too little effective support.
// For each focal cell type and each neighbour column, over the focal's cells,
//   centred = k - mean(k)                       (R's mean, refined, na.rm)
//   n_eff   = sum(centred^2) / max(centred^2)   (0 when the column is constant)
// and the column is zeroed for that focal's cells when n_eff < min_effective.
//
// Shapes: `kernel` is n * n_types column-major and is modified in place;
// `celltype` has n codes in [-1, n_types). `n_dropped` counts the pairs zeroed.
Status drop_sparse_kernel(Span<double> kernel, std::int64_t n, std::int64_t n_types,
                          Span<const int> celltype, double min_effective,
                          std::int64_t* n_dropped);

// Subtract, from every kernel column, the mean of that column within each
// (image, cell type) group, so the slopes are estimated on the within-group
// deviation only. The mean is R's mean(), as stats::ave() used.
//
// Shapes: `kernel` is n * n_types column-major, modified in place; `image` and
// `celltype` have n codes each; a cell whose code is negative in either keeps
// its own group.
Status centre_within_groups(Span<double> kernel, std::int64_t n, std::int64_t n_types,
                            Span<const int> image, Span<const int> celltype, int n_images,
                            int n_celltypes);

// R's scale() of one column: subtract colMeans() (a plain long double sum over
// n, not the refined mean) and divide by sqrt(sum(centred^2) / (n - 1)).
// `standard_deviation` (may be null) reports the divisor, which the caller uses
// to decide whether the column varies at all.
Status standardise_column(Span<double> values, std::int64_t n, double* standard_deviation);

// The anchor decision, over the per-cell-type clean profiles:
//   owner   = the type with the largest core mean for that gene (the first on a tie)
//   anchor  = owner is not this type, owner mean > owner_threshold, and this
//             type's core mean is below core_threshold of the owner's
//
// Shapes: `core_means` is n_types * n_genes column-major; `mask` is the same
// shape and receives 1 where a gene anchors that type; `n_anchor` has n_types
// entries.
Status anchor_mask(Span<const double> core_means, std::int64_t n_types, std::int64_t n_genes,
                   double owner_threshold, double core_threshold, Span<double> mask,
                   Span<int> n_anchor);

// Scale each row to its own maximum and clamp to one, leaving a row whose
// maximum is not finite or not positive at zero, and reading a missing entry as
// zero. This is the neighbour-kernel variance normalisation of the
// data-informed prior weights.
Status normalise_rows_to_max(Span<double> values, std::int64_t n_rows, std::int64_t n_cols);

}  // namespace pace

#endif  // PACE_PREPROCESS_HPP
