#include "larch/common.hpp"

Scheduler& DefaultScheduler() {
  static Scheduler scheduler{3};
  return scheduler;
}
