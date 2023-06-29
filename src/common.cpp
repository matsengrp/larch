#include "larch/common.hpp"

Scheduler& DefaultScheduler() {
  static Scheduler scheduler;
  return scheduler;
}
