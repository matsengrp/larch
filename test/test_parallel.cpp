#include "larch/parallel/parallel_common.hpp"

#include "test_common.hpp"

#include <atomic>
#include <thread>
#include <vector>

[[maybe_unused]] static void test_parallel_for() {
  const size_t outer_size = 4 * std::thread::hardware_concurrency();
  const size_t inner_size = 4 * std::thread::hardware_concurrency();

  std::atomic<size_t> counter{0};

  std::vector<size_t> outer(outer_size);
  std::iota(outer.begin(), outer.end(), 0);

  ParallelForEach(outer, [&](size_t) {
    std::vector<size_t> inner(inner_size);
    std::iota(inner.begin(), inner.end(), 0);

    ParallelForEach(inner, [&](size_t) { counter.fetch_add(1); });
  });

  TestAssert(counter.load() == outer_size * inner_size);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({test_parallel_for, "Parallel: recursive ParallelForEach", {"parallel"}});
