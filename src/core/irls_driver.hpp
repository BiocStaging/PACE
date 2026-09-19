// irls_driver.hpp -- the outer PQL iteration of the streaming solver, and the
// three full-panel passes it is built from.
//
// Every piece of per-iteration ARITHMETIC was already compiled. What stayed in
// R was the SEQUENCING: the `for (it in seq_len(n_iter))` loop itself, the
// early-stop decision, the tau bookkeeping, the convergence history and the
// progress lines. That cost a boundary crossing per pass per iteration -- with
// it a fresh Rcpp allocation of every full-panel matrix the pass returned (B,
// U, re_var, the two standard errors, alpha, num and den), each of which R then
// copied into the loop's own state -- and it kept the one part of the fit that
// decides WHICH iteration is the last in a language the core could not see.
//
// The loop below is a transcription of that R, not a redesign: the same calls
// in the same order, with every guard, floor and clamp where it was, so the fit
// is unchanged to the bit. fit_pace_mvpql_streaming() still carries the R loop
// and hands it the fused path, which is also the reference the transcription is
// gated against (tests/testthat/test-irls-driver.R).
//
// Nothing here touches R. The two elementary special functions the dispersion
// MLE and the tau prior need are passed in, exactly as dispersion_chunk() takes
// its LogDensity, and everything the loop would otherwise print or warn about
// goes out through IrlsLoopReporter, which the binding fills with Rprintf and a
// deferred R warning.
#ifndef PACE_IRLS_DRIVER_HPP
#define PACE_IRLS_DRIVER_HPP

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "core_types.hpp"
#include "sparse_product.hpp"
#include "count_stats.hpp"
#include "dispersion.hpp"
#include "gene_solve.hpp"

namespace pace {

// Where a chunk's ambient field comes from.
//
// CACHED is the historical behaviour: W %*% Y is built once, in R, and every
// pass reads a slice of it. That is the right trade while the cache is small --
// 0.03 GB on breast cancer, 0.12 GB on melanoma.
//
// It stops being the right trade at scale. On a 1.2M cell cohort at the full
// 5,001 gene panel the cache measures 5.2 GB, on top of 4.0 GB for the counts
// and another 4.0 GB for the untransposed copy the caller still holds. That
// combination pages on a 24 GB machine, and a fit that pages is not a fit: the
// observed cost went from 275 s an iteration to 8,589 s, essentially all of it
// waiting on disk.
//
// STREAMED recomputes each chunk's columns from W and Y instead. It costs a
// sparse product per pass per iteration and saves the whole cache.
//
// The two modes are BIT-IDENTICAL, verified on a real cohort across chunk
// sizes 16, 64 and 278 and at one and four threads. Two things make that true.
// The chunking is exact, because a sparse product is column-independent: a
// chunk's columns do not depend on how the panel is cut. And the accumulation
// matches CHOLMOD's because this file includes fp_no_contract.hpp -- without
// it clang contracts `scratch[row] += a * b` into an FMA at its default
// -ffp-contract=on, which is MORE accurate than a separate multiply and add
// and therefore disagrees, by about 1e-15 an entry and about 1e-6 in the
// coefficients once a fit has amplified it. That drift was mistaken for an
// ordering difference before the include was added.
struct AmbientSource {
  // Cached: a view of the whole n x G product, sliced per chunk.
  GeneBlock cached;
  // Streamed: the two operands, plus scratch reused across chunks so that the
  // per-chunk product does not allocate. The scratch is held by POINTER, not by
  // value, which is what lets ambient_block() take this by const reference and
  // still write the chunk into it: the pointers are const, what they address is
  // not. The buffers live in run_irls_loop, one set for the whole fit.
  bool streamed = false;
  CscView weights;      // W, n x n
  CscView counts;       // Y, n x G
  std::vector<int>* scratch_column_pointer = nullptr;
  std::vector<int>* scratch_row_index = nullptr;
  std::vector<double>* scratch_values = nullptr;
  CscView* scratch_view = nullptr;
};

// The ambient block for genes [first_gene, first_gene + width), from whichever
// source the caller configured. In cached mode this is a pointer offset; in
// streamed mode it computes the product into the scratch buffers and views it.
Status ambient_block(const AmbientSource& source, std::int64_t first_gene, std::int64_t width,
                     int n_threads, const InterruptCheck& interrupted, GeneBlock* out);


// One IRLS pass over every gene chunk: the working response, the float/double
// ridge decision, the per-gene solve and the lossless NaN repair.
//
// Per chunk R used to build the working response, hand z and w back, call the
// solve, take three matrices back and write them into B and U, reconverting the
// random-effect design on every one of those solve calls. Here the design is
// converted once by the caller, z and w are buffers reused across chunks, and
// only the finished coefficients are written out.
//
// `last_iter` forces the interior to double and fills the standard errors;
// otherwise `se_beta` and `se_u` are left untouched, which is how the caller's
// NA initialisation survives the interior iterations. `nan_genes` receives the
// 1-based numbers of the genes still non-finite after the double repair, in
// ascending order.
//
// Shapes: `x_fixed` is the dense n * p design, or empty when p == 1 and `x1`
// carries the single column; `solve_x_fixed` is always the dense n * p design
// the solve itself reads. `lam_diag` is q_total * n_genes; `beta_out` is
// p * n_genes; `u_out`, `re_var_out` and `se_u` are q_total * n_genes;
// `se_beta` is p * n_genes.
Status fit_pass1(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                 Span<const double> solve_x_fixed, std::int64_t p, const CscView& z,
                 const std::vector<SolveBlock>& solve_blocks, Span<const double> beta_in,
                 Span<const double> u_in, Span<const double> lam_diag, const GeneBlock& counts,
                 const AmbientSource& ambient, Span<const double> offset, Span<const double> rho,
                 Span<const double> alpha, Span<const double> sample_weight, bool nb2,
                 bool gaussian, bool seed_iteration, std::int64_t n, std::int64_t q_total,
                 std::int64_t n_genes, std::int64_t chunk_size, std::int64_t sub_genes,
                 int interior_precision, bool last_iter, Span<double> beta_out, Span<double> u_out,
                 Span<double> re_var_out, Span<double> se_beta, Span<double> se_u,
                 std::vector<int>* nan_genes, int n_threads, const InterruptCheck& interrupted);

// The post-solve rho accumulation and convergence metric over every gene chunk.
// Both linear predictors -- the current one and the previous iteration's -- are
// buffers reused across chunks rather than the pair of n x chunk R matrices the
// loop used to build and drop per chunk, and the running accumulators never
// leave. `num` and `den` are overwritten, not added to.
//
// Without `have_previous` the previous linear predictor is the solver's own
// seed, log(max(y, 0.5)) - offset, as rho_accumulate() builds it.
Status rho_pass(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                std::int64_t p, const CscView& z, Span<const double> beta_in,
                Span<const double> u_in, Span<const double> prev_beta, Span<const double> prev_u,
                bool have_previous, const GeneBlock& counts, const AmbientSource& ambient,
                Span<const double> offset, Span<const double> rho, Span<const double> alpha,
                Span<const double> mask, Span<const int> mask_index, std::int64_t n_mask_rows,
                bool nb2, std::int64_t n, std::int64_t n_genes, std::int64_t chunk_size,
                Span<double> num, Span<double> den, double* rel_delta_max, double* rel_delta_sum,
                std::int64_t* n_finite, std::int64_t* n_nonfinite, Span<double> tail_counts,
                int n_threads, const InterruptCheck& interrupted);

// The dispersion sweep over every gene chunk: the NB dispersion MLE per gene,
// or, on the Gaussian path, the per-gene residual variance. eta is a buffer
// reused across chunks rather than an n x chunk R matrix built and dropped per
// chunk. `n_noninteger` counts the genes whose counts are not whole numbers,
// which the caller warns about.
Status dispersion_pass(Span<const double> x1, bool x1_is_unit, Span<const double> x_fixed,
                       std::int64_t p, const CscView& z, Span<const double> beta_in,
                       Span<const double> u_in, const GeneBlock& counts, const AmbientSource& ambient,
                       Span<const double> offset, Span<const double> rho, bool nb2, bool gaussian,
                       bool zero_collapse, double max_cells, LogDensity density, bool fast_density,
                       std::int64_t n, std::int64_t n_genes, std::int64_t chunk_size,
                       Span<double> alpha, std::int64_t* n_noninteger, int n_threads,
                       const InterruptCheck& interrupted);

// ---------------------------------------------------------------------------
// The loop itself.
// ---------------------------------------------------------------------------

// Which prior the variance components are shrunk under.
//
// IMPORTANT: the order is the one the CALL SITE matches against -- the literal
// c("shared", "hierarchical", "adaptive", "half_cauchy") at
// engine-streaming.R:788 -- and NOT the order the `tau_shrinkage` argument is
// declared in, which is c("hierarchical", "shared", ...). Those two disagree in
// their first two entries. Reordering this enum to agree with the declaration
// would swap `shared` and `hierarchical` and silently change every fit that
// uses either. Change the enum only together with that match() vector.
enum class TauShrinkage { shared, hierarchical, adaptive, half_cauchy };

// One random-effect block as the tau bookkeeping sees it.
//
// R holds a block's variance components as a K_terms x K_groups matrix filled
// BY ROW from the block's slice of the flat length-q vector, so the flat slot
// of (term t, group g) is col_offset + t * n_groups + g while the matrix the
// hierarchical shrinkage reads is column-major. The two layouts disagree, which
// is why the hierarchical branch transposes rather than pointing at the slice,
// and why a flat slice handed straight to the shrinkage has been a bug before.
//
// `group_size` is the number of cells in each group, R's
// vapply(cells_by_grp_list[[b]], length, 0L), which the hierarchical weights read.
struct TauBlock {
  int col_offset = 0;
  int n_terms = 0;
  int n_groups = 0;
  int n_cols = 0;
  Span<const int> group_size;
};

// Everything the loop reads and never writes.
struct IrlsLoopInputs {
  Span<const double> x1;              // the single fixed column when p == 1
  bool x1_is_unit = false;            // ... and it is all ones, so the multiply is skipped
  Span<const double> x_fixed;         // the dense n * p design, empty when p == 1
  Span<const double> solve_x_fixed;   // the dense n * p design the solve reads, always present
  std::int64_t p = 0;
  CscView z;                          // the n x q random-effect design
  const std::vector<SolveBlock>* solve_blocks = nullptr;
  const std::vector<TauBlock>* tau_blocks = nullptr;
  GeneBlock counts;                   // the n x G counts, read one gene at a time
  AmbientSource ambient;              // the n x G E^tech field, cached or streamed
  Span<const double> offset;          // n entries
  Span<const double> sample_weight;   // n entries, or empty
  Span<const double> mask;            // the anchor mask, or empty
  Span<const int> mask_index;         // each cell's 1-based type, or empty for an n x G mask
  std::int64_t n_mask_rows = 0;
  Span<const double> data_informed_weights;  // q * G, or empty
  Span<const double> alpha_init;      // G entries: all ones (NB) or the sigma2 seed (Gaussian)
  std::int64_t n = 0;
  std::int64_t q = 0;
  std::int64_t n_genes = 0;
};

// Every knob the loop is given, with the same names and meanings the R
// arguments carry.
struct IrlsLoopOptions {
  std::int64_t n_iter = 16;
  std::int64_t min_iter = 12;
  double early_stop_tol = 2e-2;       // 0 disables the early stop
  double alpha_warmup = 10;           // infinite always re-fits the dispersion
  // How R renders `alpha_warmup` in the "alpha FROZEN" line. R's format() gives
  // "6" for 6 and "Inf" for Inf, and reproducing that here would be guessing at
  // R's own rules, so the caller passes the rendered text.
  std::string alpha_warmup_label;
  bool nb2 = false;
  bool gaussian = false;
  bool zero_collapse = true;
  double alpha_max_cells = 0;         // non-finite keeps every cell
  bool fast_density = true;
  int interior_precision = 1;
  std::int64_t chunk_size = 128;
  std::int64_t sub_genes = 16;
  double tau_max = 100;               // non-finite means no cap, as .clamp_tau() reads it
  TauShrinkage tau_shrinkage = TauShrinkage::adaptive;
  double d0_min = 1;                  // R_D0_MIN
  bool use_reml = false;              // R_REML_TAU
  bool rd_diag = false;               // R_RD_DIAG
  bool verbose = true;
  int n_threads = 1;
};

// Where the loop's own text goes. The core never prints and never warns; the
// binding fills these with Rprintf and a deferred R warning.
//
// `warn_nan_genes` is separate because its message names the genes, and the
// gene names are R's. It receives the iteration and the 1-based gene numbers.
struct IrlsLoopReporter {
  std::function<void(const std::string&)> message;
  std::function<void(const std::string&)> warn;
  std::function<void(std::int64_t, Span<const int>)> warn_nan_genes;
};

// Where the loop's results go. Every span is caller-owned and fully written
// except `se_beta` and `se_u`, which are written on the last iteration only and
// must arrive holding whatever the caller wants an unfinished fit to report.
//
// The history spans hold n_iter columns; only the first `*iterations_run` of
// them are written, which is what the caller returns.
struct IrlsLoopOutputs {
  Span<double> beta;                  // p * G
  Span<double> u;                     // q * G
  Span<double> se_beta;               // p * G
  Span<double> se_u;                  // q * G
  Span<double> alpha;                 // G
  Span<double> tau_g_array;           // q * G
  Span<double> tau_flat;              // q: the per-block components the fit reports
  Span<double> rho;                   // n
  Span<double> history_tau_flat;      // q * n_iter
  Span<double> history_alpha;         // G * n_iter
  Span<double> history_rel_delta;     // n_iter
  Span<double> history_rel_delta_mean;
  Span<double> history_n_nan_genes;
  Span<double> history_n_nonfinite;
  Span<double> history_tau_max_seen;
  std::int64_t* iterations_run = nullptr;
  bool* converged = nullptr;
};

// Run the outer PQL iteration.
//
// Per iteration, in this order: decide the early stop from the PREVIOUS
// iteration's mean rel_delta (the one-iteration lag is deliberate -- it makes
// the stopping iteration the last one, so it runs the standard errors in double
// and a final dispersion update, and getting it wrong changes which iteration is
// last and therefore the answer); solve the gene chunks; accumulate and shrink
// the per-cell contamination loading (negative binomial only, the Gaussian path
// taking a coefficient-based convergence metric instead); re-fit the dispersion
// unless the warm-up has frozen it; update, shrink, weight and cap the variance
// components; record the history and, when asked, print the iteration line.
//
// `density` is R's dnbinom_mu and `trigamma` is R's trigamma, passed in so the
// numbers match the optimize()/uniroot() calls they replace while the core keeps
// no R dependency of its own.
Status run_irls_loop(const IrlsLoopInputs& inputs, const IrlsLoopOptions& options,
                     LogDensity density, Trigamma trigamma, const IrlsLoopReporter& reporter,
                     const IrlsLoopOutputs& outputs, const InterruptCheck& interrupted);

}  // namespace pace

#endif  // PACE_IRLS_DRIVER_HPP
