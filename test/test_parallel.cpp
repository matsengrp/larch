#include "test_common.hpp"

#include "larch/parallel/parallel_common.hpp"
#include "larch/parallel/reduction.hpp"

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

[[maybe_unused]] static void test_parallel_reduction() {
  const size_t num_iterations = 10;
  const size_t num_items = 10 * std::thread::hardware_concurrency();
  const size_t num_buckets = 4;

  for (size_t iter = 0; iter < num_iterations; ++iter) {
    std::vector<size_t> items(num_items);
    std::iota(items.begin(), items.end(), 0);

    // First ParallelForEach - like BuildConnectionsRaw
    std::vector<std::atomic<size_t>> counters(num_items);
    ParallelForEach(items, [&](size_t item) { counters[item].store(0); });

    // Second ParallelForEach with Reduction - like BuildRootAndLeafs
    Reduction<std::vector<size_t>> reduction{num_buckets};
    ParallelForEach(items, [&](size_t item) {
      reduction.AddElement(
          [](std::vector<size_t>& vec, size_t val) { vec.push_back(val); }, item);
    });

    size_t total = 0;
    reduction.GatherAndClear([&total](auto buckets) {
      for (auto&& bucket : buckets) {
        total += bucket.size();
      }
    });

    TestAssert(total == num_items);
  }
}

[[maybe_unused]] static const auto test_added1 =
    add_test({test_parallel_reduction, "Parallel: reduction", {"parallel"}});

// Test that mimics BuildConnections pattern more closely
[[maybe_unused]] static void test_parallel_build_connections_pattern() {
  const size_t num_iterations = 50;
  const size_t num_nodes = 1000;

  for (size_t iter = 0; iter < num_iterations; ++iter) {
    std::vector<size_t> nodes(num_nodes);
    std::iota(nodes.begin(), nodes.end(), 0);

    // Simulate BuildConnectionsRaw: clear connections
    std::vector<std::vector<size_t>> connections(num_nodes);
    ParallelForEach(nodes, [&](size_t node) { connections[node].clear(); });

    // Simulate BuildRootAndLeafs: find root and collect leafs
    std::atomic<size_t> root_id{static_cast<size_t>(-1)};
    Reduction<std::vector<size_t>> leafs{32};

    ParallelForEach(nodes, [&](size_t node_id) {
      // Simulate IsUA check (node 0 is root)
      if (node_id == 0) {
        root_id.store(node_id);
      }
      // Simulate IsLeaf check (nodes > num_nodes/2 are leafs)
      if (node_id > num_nodes / 2) {
        leafs.AddElement(
            [](std::vector<size_t>& ls, size_t id) { ls.push_back(id); }, node_id);
      }
    });

    TestAssert(root_id.load() == 0);

    size_t leaf_count = 0;
    leafs.GatherAndClear([&leaf_count](auto buckets) {
      for (auto&& bucket : buckets) {
        leaf_count += bucket.size();
      }
    });

    TestAssert(leaf_count == num_nodes - num_nodes / 2 - 1);
  }
}

[[maybe_unused]] static const auto test_added2 = add_test(
    {test_parallel_build_connections_pattern, "Parallel: build connections pattern", {"parallel"}});

// Test using actual BuildConnections
[[maybe_unused]] static void test_parallel_actual_build_connections() {
  const size_t num_iterations = 20;

  for (size_t iter = 0; iter < num_iterations; ++iter) {
    auto dag = LoadDAGFromProtobuf("data/test_5_trees/tree_0.pb.gz");
    dag.View().BuildConnections();
  }
}

[[maybe_unused]] static const auto test_added3 = add_test(
    {test_parallel_actual_build_connections, "Parallel: actual BuildConnections", {"parallel"}});
