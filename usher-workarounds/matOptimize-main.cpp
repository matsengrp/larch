#include <cstdint>
#include <atomic>

std::atomic_bool interrupted(false);
bool use_bound = false;
int process_count = 0;
int this_rank = 0;
uint32_t num_threads = 0;
