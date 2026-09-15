// core_types.hpp -- shared types for the R-free PACE core.
//
// Nothing in src/core/ includes R or Rcpp headers. Functions report failure
// through Status and never throw across the core API.
#ifndef PACE_CORE_TYPES_HPP
#define PACE_CORE_TYPES_HPP

#include <cstdint>
#include <functional>
#include <string>
#include <utility>

namespace pace {

// Non-owning view of a contiguous array. The caller (usually R) owns the memory
// and keeps it alive for the duration of the call.
template <typename T>
struct Span {
  T* data = nullptr;
  std::int64_t size = 0;

  Span() = default;
  Span(T* data_in, std::int64_t size_in) : data(data_in), size(size_in) {}

  T& operator[](std::int64_t index) const { return data[index]; }
};

enum class StatusCode { ok, invalid_argument, interrupted, internal_error };

// Result of a core call: `ok`, or a code with a human-readable message.
struct Status {
  StatusCode code = StatusCode::ok;
  std::string message;

  static Status success() { return Status(); }

  static Status failure(StatusCode failure_code, std::string failure_message) {
    Status status;
    status.code = failure_code;
    status.message = std::move(failure_message);
    return status;
  }

  bool is_ok() const { return code == StatusCode::ok; }
};

// Returns true when the user has asked to stop. It is only ever called on the
// thread that started the work (the R main thread), never from a worker.
using InterruptCheck = std::function<bool()>;

}  // namespace pace

#endif  // PACE_CORE_TYPES_HPP
