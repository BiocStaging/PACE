#ifndef PACE_CORE_STAGE_TIMER_HPP
#define PACE_CORE_STAGE_TIMER_HPP

#include <atomic>
#include <chrono>

namespace pace {

// Wall-clock spent in each stage of a fit, accumulated across every call.
//
// The package reports one number per iteration and nothing below it, so every
// question about where the time goes has had to be answered by arithmetic on
// flop counts and cache-line traffic. That has been wrong often enough to be
// worth a few clock reads: a fit makes some hundreds of solve calls and some
// thousands of pass calls, so the overhead here is microseconds against hours.
//
// Counters are atomic because the passes are entered from the calling thread
// but the stages inside them run on the pool. Only whole stages are timed, from
// the calling thread, so no counter is touched by a worker.
struct StageTimings {
  std::atomic<double> solve_stage1{0.0};   // the within-block group tensors
  std::atomic<double> solve_stage2{0.0};   // the cross-block tensors
  std::atomic<double> solve_stage3{0.0};   // one Schur solve per gene
  std::atomic<double> eta{0.0};            // rebuilding the linear predictor
  std::atomic<double> working_response{0.0};
  std::atomic<double> rho{0.0};
  std::atomic<double> dispersion{0.0};
  std::atomic<double> ambient{0.0};        // the E^tech block, cached or streamed
  // parallel_for spawns and joins fresh std::threads on every call, so the
  // dispatch count is a cost in its own right, not just bookkeeping.
  std::atomic<long long> dispatches{0};
  std::atomic<double> dispatch_overhead{0.0};

  void reset() {
    solve_stage1 = 0.0; solve_stage2 = 0.0; solve_stage3 = 0.0;
    eta = 0.0; working_response = 0.0; rho = 0.0; dispersion = 0.0; ambient = 0.0;
    dispatches = 0; dispatch_overhead = 0.0;
  }
};

StageTimings& stage_timings();

// Adds its lifetime to one counter. Accumulating into an atomic<double> needs
// a compare-exchange loop: there is no atomic fetch_add for floating point
// before C++20, and these are entered from more than one thread over a fit.
class ScopedStageTimer {
 public:
  explicit ScopedStageTimer(std::atomic<double>& total)
      : total_(total), start_(std::chrono::steady_clock::now()) {}

  ~ScopedStageTimer() {
    const std::chrono::duration<double> elapsed =
        std::chrono::steady_clock::now() - start_;
    double current = total_.load(std::memory_order_relaxed);
    while (!total_.compare_exchange_weak(current, current + elapsed.count(),
                                         std::memory_order_relaxed)) {
    }
  }

  ScopedStageTimer(const ScopedStageTimer&) = delete;
  ScopedStageTimer& operator=(const ScopedStageTimer&) = delete;

 private:
  std::atomic<double>& total_;
  std::chrono::steady_clock::time_point start_;
};

}  // namespace pace

#endif  // PACE_CORE_STAGE_TIMER_HPP
