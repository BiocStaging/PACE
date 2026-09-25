// thread_pool.hpp -- deterministic parallel-for on std::thread (no OpenMP).
#ifndef PACE_THREAD_POOL_HPP
#define PACE_THREAD_POOL_HPP

#include <cstdint>
#include <functional>

#include "core_types.hpp"

namespace pace {

// Run `body(begin, end)` over the half-open item range [0, n_items), split into
// consecutive blocks of `block_size` items.
//
// Contract:
//   * Every block is processed exactly once, by one worker.
//   * `body` must write only to memory owned by the items of its block (per-item
//     output slots), so results do not depend on `n_threads`.
//   * `body` must not call R. It may allocate C++ memory.
//   * With n_threads <= 1 the blocks run in order on the calling thread.
//   * Otherwise `n_threads` workers run the blocks while the calling thread waits,
//     polling `interrupted` about every 100 ms.
//   * An exception thrown by `body` is caught inside the worker; remaining blocks
//     are skipped and Status internal_error is returned. An interrupt skips the
//     remaining blocks and returns Status interrupted.
Status parallel_for(std::int64_t n_items, int n_threads, std::int64_t block_size,
                    const std::function<void(std::int64_t, std::int64_t)>& body,
                    const InterruptCheck& interrupted);

// As above, but `body(begin, end, worker)` also receives the index of the worker
// running the block, in [0, worker_count(n_threads, n_items, block_size)).
//
// This exists so per-worker scratch can be owned by the CALLER -- one slot per
// worker, allocated before the loop -- instead of living in function-local
// `static thread_local` storage. Those destruct at every thread exit, and since
// this pool spawns and joins fresh threads on every call, a fit ran hundreds of
// thread teardowns each destroying thread_local non-POD objects inside a
// dynamically loaded library. That is a known way to lose a process silently on
// MinGW-w64, and the Windows builds died with no error, no testthat tally and no
// "Execution halted" at the first call that used more than one thread.
//
// The block partition is unchanged, so results still do not depend on n_threads.
Status parallel_for(std::int64_t n_items, int n_threads, std::int64_t block_size,
                    const std::function<void(std::int64_t, std::int64_t, int)>& body,
                    const InterruptCheck& interrupted);

// How many workers parallel_for will use for this shape. Callers size their
// per-worker scratch with this, so the two cannot disagree.
int worker_count(int n_threads, std::int64_t n_items, std::int64_t block_size);

}  // namespace pace

#endif  // PACE_THREAD_POOL_HPP
