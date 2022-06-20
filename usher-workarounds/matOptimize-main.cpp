#include <cstdint>
#include <atomic>

std::atomic_bool interrupted(false);
bool use_bound;
int process_count;
int this_rank;
uint32_t num_threads;
