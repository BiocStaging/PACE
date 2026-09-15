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

}  // namespace pace

#endif  // PACE_THREAD_POOL_HPP
