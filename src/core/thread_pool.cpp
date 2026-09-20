// thread_pool.cpp -- implementation of pace::parallel_for (see thread_pool.hpp).
#include "thread_pool.hpp"
#include "stage_timer.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace pace {

namespace {

// Number of blocks needed to cover n_items with blocks of block_size.
std::int64_t count_blocks(std::int64_t n_items, std::int64_t block_size) {
  return (n_items + block_size - 1) / block_size;
}

}  // namespace

Status parallel_for(std::int64_t n_items, int n_threads, std::int64_t block_size,
                    const std::function<void(std::int64_t, std::int64_t)>& body,
                    const InterruptCheck& interrupted) {
  if (n_items <= 0) return Status::success();
  if (block_size < 1) block_size = 1;
  const std::int64_t n_blocks = count_blocks(n_items, block_size);

  // Serial path: blocks in order on the calling thread, interrupt checks between blocks.
  if (n_threads <= 1 || n_blocks == 1) {
    for (std::int64_t block = 0; block < n_blocks; ++block) {
      if (interrupted && interrupted()) {
        return Status::failure(StatusCode::interrupted, "interrupted by the user");
      }
      const std::int64_t begin = block * block_size;
      const std::int64_t end = std::min(n_items, begin + block_size);
      try {
        body(begin, end);
      } catch (const std::exception& error) {
        return Status::failure(StatusCode::internal_error, error.what());
      } catch (...) {
        return Status::failure(StatusCode::internal_error, "unknown error in worker");
      }
    }
    return Status::success();
  }

  // Never more workers than blocks, nor (when known) than hardware threads. The
  // block partition does not depend on the worker count, so neither do results.
  std::int64_t worker_cap = std::min<std::int64_t>(n_threads, n_blocks);
  const unsigned hardware_threads = std::thread::hardware_concurrency();
  if (hardware_threads > 0) worker_cap = std::min<std::int64_t>(worker_cap, hardware_threads);
  const int n_workers = static_cast<int>(std::max<std::int64_t>(worker_cap, 1));
  std::atomic<std::int64_t> next_block{0};
  std::atomic<bool> stop{false};
  std::mutex state_mutex;
  std::condition_variable workers_done;
  int finished_workers = 0;
  std::string error_message;

  auto worker = [&]() {
    while (!stop.load()) {
      const std::int64_t block = next_block.fetch_add(1);
      if (block >= n_blocks) break;
      const std::int64_t begin = block * block_size;
      const std::int64_t end = std::min(n_items, begin + block_size);
      try {
        body(begin, end);
      } catch (const std::exception& error) {
        std::lock_guard<std::mutex> lock(state_mutex);
        if (error_message.empty()) error_message = error.what();
        stop.store(true);
      } catch (...) {
        std::lock_guard<std::mutex> lock(state_mutex);
        if (error_message.empty()) error_message = "unknown error in worker";
        stop.store(true);
      }
    }
    std::lock_guard<std::mutex> lock(state_mutex);
    finished_workers += 1;
    workers_done.notify_all();
  };

  // Spawn and join, timed: this is what a persistent pool would remove.
  stage_timings().dispatches.fetch_add(1, std::memory_order_relaxed);
  const std::chrono::steady_clock::time_point dispatch_start =
      std::chrono::steady_clock::now();

  std::vector<std::thread> threads;
  threads.reserve(n_workers);
  try {
    for (int t = 0; t < n_workers; ++t) threads.emplace_back(worker);
  } catch (const std::exception& error) {
    // Could not start every worker: stop and join the ones that did start, so no
    // joinable std::thread is destroyed (which would call std::terminate).
    stop.store(true);
    for (auto& thread : threads) thread.join();
    return Status::failure(StatusCode::internal_error,
                           std::string("could not start worker threads: ") + error.what());
  }

  bool was_interrupted = false;
  {
    std::unique_lock<std::mutex> lock(state_mutex);
    while (finished_workers < n_workers) {
      workers_done.wait_for(lock, std::chrono::milliseconds(100));
      if (finished_workers >= n_workers) break;
      if (!was_interrupted && interrupted) {
        lock.unlock();
        const bool stop_requested = interrupted();
        lock.lock();
        if (stop_requested) {
          was_interrupted = true;
          stop.store(true);
        }
      }
    }
  }
  for (auto& thread : threads) thread.join();
  {
    const std::chrono::duration<double> spent =
        std::chrono::steady_clock::now() - dispatch_start;
    std::atomic<double>& total = stage_timings().dispatch_overhead;
    double current = total.load(std::memory_order_relaxed);
    while (!total.compare_exchange_weak(current, current + spent.count(),
                                        std::memory_order_relaxed)) {
    }
  }

  if (!error_message.empty()) return Status::failure(StatusCode::internal_error, error_message);
  if (was_interrupted) return Status::failure(StatusCode::interrupted, "interrupted by the user");
  return Status::success();
}

}  // namespace pace
