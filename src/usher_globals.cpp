#ifdef USE_USHER

#include "larch/usher_glue.hpp"
#include <tbb/task_scheduler_init.h>

std::atomic_bool interrupted(false);
int process_count = 1;
int this_rank = 0;
uint32_t num_threads =
    static_cast<uint32_t>(tbb::task_scheduler_init::default_num_threads());
FILE* movalbe_src_log;
bool changing_radius = false;
bool use_bound = true;

#endif
