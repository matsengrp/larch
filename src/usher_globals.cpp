#ifdef USE_USHER

#include "larch/usher_glue.hpp"
#include <tbb/info.h>

std::atomic_bool interrupted(false);
int process_count = 1;
int this_rank = 0;
uint32_t num_threads =
    static_cast<uint32_t>(tbb::info::default_concurrency());
FILE* movalbe_src_log;
bool changing_radius = false;
bool use_bound = true;

#endif
