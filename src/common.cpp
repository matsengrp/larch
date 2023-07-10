#include "larch/parallel/for_loop.hpp"

Scheduler& DefaultScheduler() {
  static Scheduler scheduler;
  return scheduler;
}
