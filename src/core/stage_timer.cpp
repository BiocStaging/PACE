#include "stage_timer.hpp"

namespace pace {

StageTimings& stage_timings() {
  static StageTimings timings;
  return timings;
}

}  // namespace pace
