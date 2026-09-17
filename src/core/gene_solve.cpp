// gene_solve.cpp -- implementation of gene_solve.hpp.
//
// The arithmetic is carried over unchanged from src/solve_chunk_full.cpp, which
// this replaces: same three stages, same expressions, same Eigen calls, so the
// fits it produces are the fits that file produced. What changed is the shape of
// the code, not its numbers: no R types cross into it, and the two parallel
// loops run on the deterministic std::thread pool instead of OpenMP.
#include "gene_solve.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <memory>
#include <limits>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <vector>

#include "thread_pool.hpp"

namespace pace {
namespace {

// Eigen opens a thread team of its own for a large product unless it is already
// inside one. The kernel this replaces called stages 1 and 3 from inside an
// OpenMP parallel region, where Eigen serialises itself, and stage 2 from the
// calling thread, where it does not. Here the callers are std::threads, which
// Eigen cannot see, so the same division is made explicitly: without it the
// nested teams oversubscribe the machine (a whole fit took five times as long)
// and the products come out differently from the kernel this has to reproduce.
class EigenThreads {
 public:
  explicit EigenThreads(int threads) : previous_(Eigen::nbThreads()) {
    Eigen::setNbThreads(threads);
  }
  ~EigenThreads() { Eigen::setNbThreads(previous_); }
  EigenThreads(const EigenThreads&) = delete;
  EigenThreads& operator=(const EigenThreads&) = delete;

 private:
  int previous_;
};

// The cross-block tensor of one (group in block 1, group in block 2) pair.
template <typename T>
struct CrossBlock {
  int n_terms_1 = 0;
  int n_terms_2 = 0;
  int n_groups_1 = 0;
  int n_groups_2 = 0;
  int group_1 = 0;
  int group_2 = 0;
  int col_offset_1 = 0;
  int col_offset_2 = 0;
  Eigen::Matrix<T, -1, -1> tensor;  // (n_terms_1 * n_terms_2) x n_genes
};

template <typename T>
Eigen::Matrix<T, -1, -1> to_matrix(Span<const double> values, std::int64_t rows,
                                   std::int64_t cols) {
  Eigen::Matrix<T, -1, -1> out(rows, cols);
  for (std::int64_t j = 0; j < cols; ++j) {
    for (std::int64_t i = 0; i < rows; ++i) {
      out(i, j) = static_cast<T>(values[i + j * rows]);
    }
  }
  return out;
}

template <typename T>
Status solve_impl(Span<const double> x_fixed_in, std::int64_t n, std::int64_t p,
                  const std::vector<SolveBlock>& blocks, Span<const double> w_in,
                  Span<const double> z_in, Span<const double> lam_in, std::int64_t q_total,
                  std::int64_t n_genes, int n_threads, const InterruptCheck& interrupted,
                  Span<double> beta_out, Span<double> u_out, Span<double> ainv_diag_out) {
  using Matrix = Eigen::Matrix<T, -1, -1>;
  using Vector = Eigen::Matrix<T, -1, 1>;
  const int n_blocks = static_cast<int>(blocks.size());

  const Matrix x_fixed = to_matrix<T>(x_fixed_in, n, p);
  const Matrix w = to_matrix<T>(w_in, n, n_genes);
  const Matrix z = to_matrix<T>(z_in, n, n_genes);
  const Matrix lam = to_matrix<T>(lam_in, q_total, n_genes);
  std::vector<Matrix> terms(n_blocks);
  for (int b = 0; b < n_blocks; ++b) {
    terms[b] = to_matrix<T>(blocks[b].terms, n, blocks[b].n_terms);
  }

  // ---- Stage 1: the within-block tensors of each (block, group) ----
  // The core's own workers carry the parallelism from here on.
  std::unique_ptr<EigenThreads> serial_eigen(new EigenThreads(1));
  std::vector<std::vector<Matrix>> ztwz(n_blocks);
  std::vector<std::vector<Matrix>> xtwz(n_blocks);
  std::vector<std::vector<Matrix>> ztwz_rhs(n_blocks);
  for (int b = 0; b < n_blocks; ++b) {
    const int n_terms = blocks[b].n_terms;
    const int n_groups = blocks[b].n_groups;
    ztwz[b].resize(n_groups);
    xtwz[b].resize(n_groups);
    ztwz_rhs[b].resize(n_groups);
    const Matrix& block_terms = terms[b];
    // Groups run one at a time and the GENES are what is parallelised.
    //
    // Parallelising over groups instead cannot balance: a group's work is
    // proportional to its cell count, and one cell type routinely holds most of
    // the cohort (81.5% in the melanoma cohort, 69.6% in the 1.2M melanoma TMA).
    // The makespan is then that single group no matter how many workers there
    // are, which is why this stage measured 19.9 s on four threads and 18.5 s on
    // eight while the dispersion step in the same run scaled 1.51x.
    //
    // Tiles are multiples of Eigen's gebp register width (Traits::nr == 4 for
    // float and double on this target) and any ragged tail is kept as its own
    // final tile, so every column goes through the same nr-wide kernel path it
    // would have taken in one big multiply. Eigen's serial blocking picks kc and
    // mc independently of the number of columns, and its nc loop sits inside the
    // k loop, so an output column sees the same sequence of partial sums either
    // way: the result is bit-identical, which a test asserts.
    for (int g = 0; g < n_groups; ++g) {
        const int first = blocks[b].group_offsets[g];
        const int last = blocks[b].group_offsets[g + 1];
        const int n_g = last - first;
        if (n_g == 0) {
          ztwz[b][g] = Matrix::Zero(n_terms * n_terms, n_genes);
          xtwz[b][g] = Matrix::Zero(p * n_terms, n_genes);
          ztwz_rhs[b][g] = Matrix::Zero(n_terms, n_genes);
          continue;
        }
        Matrix terms_g(n_g, n_terms);
        Matrix x_g(n_g, p);
        for (int i = 0; i < n_g; ++i) {
          const int row = blocks[b].group_cells[first + i];
          terms_g.row(i) = block_terms.row(row);
          x_g.row(i) = x_fixed.row(row);
        }
        Matrix term_pairs(n_g, n_terms * n_terms);
        for (int t2 = 0; t2 < n_terms; ++t2) {
          for (int t1 = 0; t1 < n_terms; ++t1) {
            term_pairs.col(t2 * n_terms + t1) = terms_g.col(t1).cwiseProduct(terms_g.col(t2));
          }
        }
        Matrix fixed_term_pairs(n_g, p * n_terms);
        for (int t = 0; t < n_terms; ++t) {
          for (std::int64_t pi = 0; pi < p; ++pi) {
            fixed_term_pairs.col(t * p + pi) = x_g.col(pi).cwiseProduct(terms_g.col(t));
          }
        }
        ztwz[b][g].resize(n_terms * n_terms, n_genes);
        xtwz[b][g].resize(p * n_terms, n_genes);
        ztwz_rhs[b][g].resize(n_terms, n_genes);

        // Only large groups are tiled. Eigen's serial blocking picks mc
        // independently of the column count ONLY once k (here n_g) is big
        // enough to be k-blocked; below that mc moves with the column count and
        // tiling would change the arithmetic. Small groups also carry too
        // little work for the threading to pay. The fixtures are 2k-9k cells
        // over 6-13 types, so their groups sit in the low hundreds and take the
        // untiled path -- which is why those configurations stay bit-identical.
        const int min_cells_to_tile = 4096;
        if (n_g < min_cells_to_tile) {
          Matrix w_g(n_g, n_genes);
          Matrix z_g(n_g, n_genes);
          for (std::int64_t j = 0; j < n_genes; ++j) {
            const T* w_source = w.col(j).data();
            const T* z_source = z.col(j).data();
            T* w_target = w_g.col(j).data();
            T* z_target = z_g.col(j).data();
            for (int i = 0; i < n_g; ++i) {
              const int row = blocks[b].group_cells[first + i];
              w_target[i] = w_source[row];
              z_target[i] = z_source[row];
            }
          }
          ztwz[b][g].noalias() = term_pairs.transpose() * w_g;
          xtwz[b][g].noalias() = fixed_term_pairs.transpose() * w_g;
          const Matrix weighted_z = w_g.cwiseProduct(z_g);
          ztwz_rhs[b][g].noalias() = terms_g.transpose() * weighted_z;
          continue;
        }

        // Column tiles: multiples of four, with the ragged tail left whole.
        const std::int64_t register_width = 4;
        const std::int64_t aligned = (n_genes / register_width) * register_width;
        std::vector<std::pair<std::int64_t, std::int64_t>> tiles;
        const std::int64_t tile_width = register_width * 4;
        for (std::int64_t start = 0; start < aligned; start += tile_width) {
          tiles.emplace_back(start, std::min(tile_width, aligned - start));
        }
        if (aligned < n_genes) tiles.emplace_back(aligned, n_genes - aligned);

        const int* group_cells = blocks[b].group_cells.data + first;
        auto tile_body = [&](std::int64_t begin, std::int64_t end) {
          for (std::int64_t t = begin; t < end; ++t) {
            const std::int64_t column = tiles[static_cast<std::size_t>(t)].first;
            const std::int64_t width = tiles[static_cast<std::size_t>(t)].second;
            Matrix w_tile(n_g, width);
            Matrix z_tile(n_g, width);
            for (std::int64_t j = 0; j < width; ++j) {
              const T* w_source = w.col(column + j).data();
              const T* z_source = z.col(column + j).data();
              T* w_target = w_tile.col(j).data();
              T* z_target = z_tile.col(j).data();
              for (int i = 0; i < n_g; ++i) {
                w_target[i] = w_source[group_cells[i]];
                z_target[i] = z_source[group_cells[i]];
              }
            }
            ztwz[b][g].middleCols(column, width).noalias() =
                term_pairs.transpose() * w_tile;
            xtwz[b][g].middleCols(column, width).noalias() =
                fixed_term_pairs.transpose() * w_tile;
            const Matrix weighted_z = w_tile.cwiseProduct(z_tile);
            ztwz_rhs[b][g].middleCols(column, width).noalias() =
                terms_g.transpose() * weighted_z;
          }
        };
        const Status tile_status = parallel_for(static_cast<std::int64_t>(tiles.size()),
                                                n_threads, 1, tile_body, interrupted);
        if (!tile_status.is_ok()) return tile_status;
    }
  }

  // ---- Stage 2: the cross-block tensors, over the cells two groups share ----
  // On the calling thread, with Eigen's own parallelism, as before.
  serial_eigen.reset();
  std::vector<CrossBlock<T>> cross;
  for (int b1 = 0; b1 + 1 < n_blocks; ++b1) {
    const int n_terms_1 = blocks[b1].n_terms;
    for (int b2 = b1 + 1; b2 < n_blocks; ++b2) {
      const int n_terms_2 = blocks[b2].n_terms;
      for (int g1 = 0; g1 < blocks[b1].n_groups; ++g1) {
        const int first = blocks[b1].group_offsets[g1];
        const int last = blocks[b1].group_offsets[g1 + 1];
        if (last == first) continue;
        std::unordered_map<int, std::vector<int>> buckets;
        for (int k = first; k < last; ++k) {
          const int cell = blocks[b1].group_cells[k];
          buckets[blocks[b2].cell_group[cell] - 1].push_back(cell);
        }
        for (const auto& bucket : buckets) {
          const std::vector<int>& cells = bucket.second;
          const int n_shared = static_cast<int>(cells.size());
          if (n_shared == 0) continue;
          Matrix terms_1(n_shared, n_terms_1);
          Matrix terms_2(n_shared, n_terms_2);
          Matrix w_shared(n_shared, n_genes);
          for (int i = 0; i < n_shared; ++i) {
            terms_1.row(i) = terms[b1].row(cells[i]);
            terms_2.row(i) = terms[b2].row(cells[i]);
            w_shared.row(i) = w.row(cells[i]);
          }
          Matrix products(n_shared, n_terms_1 * n_terms_2);
          for (int t2 = 0; t2 < n_terms_2; ++t2) {
            for (int t1 = 0; t1 < n_terms_1; ++t1) {
              products.col(t2 * n_terms_1 + t1) = terms_1.col(t1).cwiseProduct(terms_2.col(t2));
            }
          }
          CrossBlock<T> entry;
          entry.n_terms_1 = n_terms_1;
          entry.n_terms_2 = n_terms_2;
          entry.n_groups_1 = blocks[b1].n_groups;
          entry.n_groups_2 = blocks[b2].n_groups;
          entry.group_1 = g1;
          entry.group_2 = bucket.first;
          entry.col_offset_1 = blocks[b1].col_offset;
          entry.col_offset_2 = blocks[b2].col_offset;
          entry.tensor.noalias() = products.transpose() * w_shared;
          cross.push_back(std::move(entry));
        }
      }
    }
  }

  // ---- Stage 3: one solve per gene ----
  serial_eigen.reset(new EigenThreads(1));
  const T missing = std::numeric_limits<T>::quiet_NaN();
  for (std::int64_t k = 0; k < p * n_genes; ++k) beta_out[k] = missing;
  for (std::int64_t k = 0; k < q_total * n_genes; ++k) u_out[k] = missing;
  for (std::int64_t k = 0; k < (p + q_total) * n_genes; ++k) ainv_diag_out[k] = missing;

  const bool use_schur = (n_blocks == 2);
  auto body = [&](std::int64_t begin, std::int64_t end) {
    for (std::int64_t gi = begin; gi < end; ++gi) {
      const Vector w_gene = w.col(gi);
      const Vector z_gene = z.col(gi);
      const Matrix weighted_x = x_fixed.array().colwise() * w_gene.array();
      const Matrix xtwx = x_fixed.transpose() * weighted_x;
      const Vector xtwz_rhs = x_fixed.transpose() * (w_gene.cwiseProduct(z_gene));

      if (use_schur) {
        const int n_terms_1 = blocks[0].n_terms;
        const int n_groups_1 = blocks[0].n_groups;
        const int n_terms_2 = blocks[1].n_terms;
        const int n_groups_2 = blocks[1].n_groups;
        const std::int64_t q1 = static_cast<std::int64_t>(n_terms_1) * n_groups_1;
        const std::int64_t q2 = static_cast<std::int64_t>(n_terms_2) * n_groups_2;
        const int offset_1 = blocks[0].col_offset;
        const int offset_2 = blocks[1].col_offset;

        std::vector<Eigen::LLT<Matrix>> first_block_chol;
        first_block_chol.reserve(n_groups_1);
        bool ok = true;
        for (int g = 0; g < n_groups_1; ++g) {
          Matrix gram(n_terms_1, n_terms_1);
          const Matrix& tensor = ztwz[0][g];
          for (int t2 = 0; t2 < n_terms_1; ++t2) {
            for (int t1 = 0; t1 < n_terms_1; ++t1) {
              gram(t1, t2) = tensor(t1 + t2 * n_terms_1, gi);
            }
          }
          for (int t = 0; t < n_terms_1; ++t) {
            gram(t, t) += lam(offset_1 + t * n_groups_1 + g, gi);
          }
          first_block_chol.emplace_back(gram);
          if (first_block_chol.back().info() != Eigen::Success) {
            ok = false;
            break;
          }
        }
        if (!ok) continue;

        std::vector<Matrix> second_block(n_groups_2);
        for (int g = 0; g < n_groups_2; ++g) {
          Matrix gram(n_terms_2, n_terms_2);
          const Matrix& tensor = ztwz[1][g];
          for (int t2 = 0; t2 < n_terms_2; ++t2) {
            for (int t1 = 0; t1 < n_terms_2; ++t1) {
              gram(t1, t2) = tensor(t1 + t2 * n_terms_2, gi);
            }
          }
          for (int t = 0; t < n_terms_2; ++t) {
            gram(t, t) += lam(offset_2 + t * n_groups_2 + g, gi);
          }
          second_block[g] = std::move(gram);
        }

        Matrix coupling = Matrix::Zero(q1, q2);
        for (const CrossBlock<T>& entry : cross) {
          const int row_offset = entry.group_1 * n_terms_1;
          const int col_offset = entry.group_2 * n_terms_2;
          for (int t2 = 0; t2 < n_terms_2; ++t2) {
            for (int t1 = 0; t1 < n_terms_1; ++t1) {
              coupling(row_offset + t1, col_offset + t2) = entry.tensor(t1 + t2 * n_terms_1, gi);
            }
          }
        }

        Matrix solved_coupling(q1, q2);
        for (int g = 0; g < n_groups_1; ++g) {
          const int offset = g * n_terms_1;
          solved_coupling.middleRows(offset, n_terms_1) =
              first_block_chol[g].solve(coupling.middleRows(offset, n_terms_1));
        }

        Matrix schur = -(coupling.transpose() * solved_coupling);
        for (int g = 0; g < n_groups_2; ++g) {
          schur.block(g * n_terms_2, g * n_terms_2, n_terms_2, n_terms_2) += second_block[g];
        }
        Eigen::LLT<Matrix> schur_chol(schur);
        if (schur_chol.info() != Eigen::Success) continue;

        Matrix xtwz_permuted(p, q_total);
        for (int g = 0; g < n_groups_1; ++g) {
          const Matrix& tensor = xtwz[0][g];
          const int col_offset = g * n_terms_1;
          for (int t = 0; t < n_terms_1; ++t) {
            for (std::int64_t pi = 0; pi < p; ++pi) {
              xtwz_permuted(pi, col_offset + t) = tensor(pi + t * p, gi);
            }
          }
        }
        for (int g = 0; g < n_groups_2; ++g) {
          const Matrix& tensor = xtwz[1][g];
          const std::int64_t col_offset = q1 + g * n_terms_2;
          for (int t = 0; t < n_terms_2; ++t) {
            for (std::int64_t pi = 0; pi < p; ++pi) {
              xtwz_permuted(pi, col_offset + t) = tensor(pi + t * p, gi);
            }
          }
        }

        Vector rhs_permuted(q_total);
        for (int g = 0; g < n_groups_1; ++g) {
          const Matrix& tensor = ztwz_rhs[0][g];
          for (int t = 0; t < n_terms_1; ++t) rhs_permuted(g * n_terms_1 + t) = tensor(t, gi);
        }
        for (int g = 0; g < n_groups_2; ++g) {
          const Matrix& tensor = ztwz_rhs[1][g];
          for (int t = 0; t < n_terms_2; ++t) rhs_permuted(q1 + g * n_terms_2 + t) = tensor(t, gi);
        }

        Vector first_solve(q1);
        for (int g = 0; g < n_groups_1; ++g) {
          const int offset = g * n_terms_1;
          first_solve.segment(offset, n_terms_1) =
              first_block_chol[g].solve(rhs_permuted.segment(offset, n_terms_1));
        }
        const Vector second_rhs = rhs_permuted.tail(q2) - coupling.transpose() * first_solve;
        const Vector effect_2 = schur_chol.solve(second_rhs);
        const Vector effect_1 = first_solve - solved_coupling * effect_2;

        Matrix first_solve_fixed(q1, p);
        for (int g = 0; g < n_groups_1; ++g) {
          const int offset = g * n_terms_1;
          first_solve_fixed.middleRows(offset, n_terms_1) =
              first_block_chol[g].solve(xtwz_permuted.block(0, offset, p, n_terms_1).transpose());
        }
        const Matrix second_rhs_fixed =
            xtwz_permuted.rightCols(q2).transpose() - coupling.transpose() * first_solve_fixed;
        const Matrix fixed_2 = schur_chol.solve(second_rhs_fixed);
        const Matrix fixed_1 = first_solve_fixed - solved_coupling * fixed_2;

        Matrix fixed_gram = xtwx;
        fixed_gram.noalias() -= xtwz_permuted.leftCols(q1) * fixed_1;
        fixed_gram.noalias() -= xtwz_permuted.rightCols(q2) * fixed_2;

        Vector fixed_rhs = xtwz_rhs;
        fixed_rhs.noalias() -= xtwz_permuted.leftCols(q1) * effect_1;
        fixed_rhs.noalias() -= xtwz_permuted.rightCols(q2) * effect_2;

        Eigen::LLT<Matrix> fixed_chol(fixed_gram);
        if (fixed_chol.info() != Eigen::Success) continue;
        const Vector beta = fixed_chol.solve(fixed_rhs);
        const Vector u_1 = effect_1 - fixed_1 * beta;
        const Vector u_2 = effect_2 - fixed_2 * beta;

        for (std::int64_t pi = 0; pi < p; ++pi) beta_out[pi + gi * p] = beta(pi);
        for (int g = 0; g < n_groups_1; ++g) {
          for (int t = 0; t < n_terms_1; ++t) {
            u_out[(offset_1 + t * n_groups_1 + g) + gi * q_total] = u_1(g * n_terms_1 + t);
          }
        }
        for (int g = 0; g < n_groups_2; ++g) {
          for (int t = 0; t < n_terms_2; ++t) {
            u_out[(offset_2 + t * n_groups_2 + g) + gi * q_total] = u_2(g * n_terms_2 + t);
          }
        }

        const Matrix fixed_inverse = fixed_chol.solve(Matrix::Identity(p, p));
        for (std::int64_t pi = 0; pi < p; ++pi) {
          ainv_diag_out[pi + gi * (p + q_total)] = fixed_inverse(pi, pi);
        }
        const Matrix schur_solved_coupling = schur_chol.solve(solved_coupling.transpose());
        Vector first_inverse_diag(q1);
        for (int g = 0; g < n_groups_1; ++g) {
          const Matrix lower_inverse =
              first_block_chol[g].matrixL().solve(Matrix::Identity(n_terms_1, n_terms_1));
          const int offset = g * n_terms_1;
          for (int t = 0; t < n_terms_1; ++t) {
            first_inverse_diag(offset + t) = lower_inverse.col(t).squaredNorm();
          }
        }
        const Matrix schur_lower_inverse =
            schur_chol.matrixL().solve(Matrix::Identity(q2, q2));
        for (int g = 0; g < n_groups_1; ++g) {
          const int offset = g * n_terms_1;
          for (int t = 0; t < n_terms_1; ++t) {
            const std::int64_t i = offset + t;
            const T random_part =
                first_inverse_diag(i) + solved_coupling.row(i).dot(schur_solved_coupling.col(i));
            const T fixed_part = fixed_1.row(i) * fixed_inverse * fixed_1.row(i).transpose();
            const std::int64_t column = offset_1 + t * n_groups_1 + g;
            ainv_diag_out[(p + column) + gi * (p + q_total)] = random_part + fixed_part;
          }
        }
        for (int g = 0; g < n_groups_2; ++g) {
          const int offset = g * n_terms_2;
          for (int t = 0; t < n_terms_2; ++t) {
            const std::int64_t i = offset + t;
            const T random_part = schur_lower_inverse.col(i).squaredNorm();
            const T fixed_part = fixed_2.row(i) * fixed_inverse * fixed_2.row(i).transpose();
            const std::int64_t column = offset_2 + t * n_groups_2 + g;
            ainv_diag_out[(p + column) + gi * (p + q_total)] = random_part + fixed_part;
          }
        }
        continue;
      }

      // ---- the dense path, for any number of blocks ----
      Matrix random_gram = Matrix::Zero(q_total, q_total);
      Matrix fixed_random = Matrix::Zero(p, q_total);
      Vector random_rhs = Vector::Zero(q_total);
      for (int b = 0; b < n_blocks; ++b) {
        const int n_terms = blocks[b].n_terms;
        const int n_groups = blocks[b].n_groups;
        const int offset = blocks[b].col_offset;
        for (int g = 0; g < n_groups; ++g) {
          const Matrix& gram_tensor = ztwz[b][g];
          const Matrix& fixed_tensor = xtwz[b][g];
          const Matrix& rhs_tensor = ztwz_rhs[b][g];
          for (int t1 = 0; t1 < n_terms; ++t1) {
            const std::int64_t col1 = offset + t1 * n_groups + g;
            for (int t2 = 0; t2 < n_terms; ++t2) {
              const std::int64_t col2 = offset + t2 * n_groups + g;
              random_gram(col1, col2) = gram_tensor(t1 + t2 * n_terms, gi);
            }
            for (std::int64_t pi = 0; pi < p; ++pi) {
              fixed_random(pi, col1) = fixed_tensor(pi + t1 * p, gi);
            }
            random_rhs(col1) = rhs_tensor(t1, gi);
          }
        }
      }
      for (const CrossBlock<T>& entry : cross) {
        for (int t1 = 0; t1 < entry.n_terms_1; ++t1) {
          const std::int64_t col1 = entry.col_offset_1 + t1 * entry.n_groups_1 + entry.group_1;
          for (int t2 = 0; t2 < entry.n_terms_2; ++t2) {
            const std::int64_t col2 = entry.col_offset_2 + t2 * entry.n_groups_2 + entry.group_2;
            const T value = entry.tensor(t1 + t2 * entry.n_terms_1, gi);
            random_gram(col1, col2) = value;
            random_gram(col2, col1) = value;
          }
        }
      }
      for (std::int64_t i = 0; i < q_total; ++i) random_gram(i, i) += lam(i, gi);

      const std::int64_t size = p + q_total;
      Matrix system(size, size);
      system.topLeftCorner(p, p) = xtwx;
      system.topRightCorner(p, q_total) = fixed_random;
      system.bottomLeftCorner(q_total, p) = fixed_random.transpose();
      system.bottomRightCorner(q_total, q_total) = random_gram;
      Vector rhs(size);
      rhs.head(p) = xtwz_rhs;
      rhs.tail(q_total) = random_rhs;

      Eigen::LLT<Matrix> chol(system);
      if (chol.info() != Eigen::Success) continue;
      const Vector solution = chol.solve(rhs);
      for (std::int64_t pi = 0; pi < p; ++pi) beta_out[pi + gi * p] = solution(pi);
      for (std::int64_t k = 0; k < q_total; ++k) u_out[k + gi * q_total] = solution(p + k);
      const Matrix lower_inverse = chol.matrixL().solve(Matrix::Identity(size, size));
      for (std::int64_t k = 0; k < size; ++k) {
        ainv_diag_out[k + gi * size] = lower_inverse.col(k).squaredNorm();
      }
    }
  };
  return parallel_for(n_genes, n_threads, 1, body, interrupted);
}

}  // namespace

Status solve_genes_chunk(Span<const double> x_fixed, std::int64_t n, std::int64_t p,
                         const std::vector<SolveBlock>& blocks, Span<const double> w,
                         Span<const double> z, Span<const double> lam_diag,
                         std::int64_t q_total, std::int64_t n_genes, bool single_precision,
                         int n_threads, const InterruptCheck& interrupted, Span<double> beta_out,
                         Span<double> u_out, Span<double> ainv_diag_out) {
  if (x_fixed.size != n * p || w.size != n * n_genes || z.size != n * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "x_fixed, w or z have the wrong size");
  }
  if (lam_diag.size != q_total * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "lam_diag must be q_total x n_genes");
  }
  if (beta_out.size != p * n_genes || u_out.size != q_total * n_genes ||
      ainv_diag_out.size != (p + q_total) * n_genes) {
    return Status::failure(StatusCode::invalid_argument, "outputs have the wrong size");
  }
  for (const SolveBlock& block : blocks) {
    if (block.terms.size != n * block.n_terms || block.cell_group.size != n ||
        block.group_offsets.size != block.n_groups + 1) {
      return Status::failure(StatusCode::invalid_argument, "a random-effect block is inconsistent");
    }
    for (std::int64_t k = 0; k < block.group_cells.size; ++k) {
      if (block.group_cells[k] < 0 || block.group_cells[k] >= n) {
        return Status::failure(StatusCode::invalid_argument, "a group cell index is out of range");
      }
    }
    for (std::int64_t i = 0; i < n; ++i) {
      if (block.cell_group[i] < 1 || block.cell_group[i] > block.n_groups) {
        return Status::failure(StatusCode::invalid_argument, "a cell group index is out of range");
      }
    }
  }
  if (single_precision) {
    return solve_impl<float>(x_fixed, n, p, blocks, w, z, lam_diag, q_total, n_genes, n_threads,
                             interrupted, beta_out, u_out, ainv_diag_out);
  }
  return solve_impl<double>(x_fixed, n, p, blocks, w, z, lam_diag, q_total, n_genes, n_threads,
                            interrupted, beta_out, u_out, ainv_diag_out);
}

}  // namespace pace
