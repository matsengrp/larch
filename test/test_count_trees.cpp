#include "subtree_weight.hpp"
#include "tree_count.hpp"

#include <iostream>
#include <string_view>

#include "test_common.hpp"

#include "dag_loader.hpp"

static void test_tree_count(MADAG& dag, size_t expected_score) {
  if (dag.GetEdgeMutations().empty()) {
    dag.GetEdgeMutations() = dag.ComputeEdgeMutations(dag.GetReferenceSequence());
  }

  SubtreeWeight<TreeCount> treecount(dag);

  size_t score = treecount.ComputeWeightBelow(dag.GetDAG().GetRoot(), {});

  assert_equal(score, expected_score,
          "True tree count " + std::to_string(expected_score) +
          " but counted " + std::to_string(score));
}

static void test_tree_count(std::string_view path, size_t expected_score) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_tree_count(dag, expected_score);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_tree_count("data/testcase/full_dag.pb.gz", 818); },
              "Subtree weight: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_tree_count("data/testcase1/full_dag.pb.gz", 7); },
              "Subtree weight: testcase1"});
