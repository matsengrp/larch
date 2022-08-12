#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>
#include <string_view>

#include "test_common.hpp"

#include "dag_loader.hpp"

static void test_subtree_weight(MADAG& dag, size_t expected_score) {
  if (dag.GetEdgeMutations().empty()) {
    dag.GetEdgeMutations() = dag.ComputeEdgeMutations(dag.GetReferenceSequence());
  }

  SubtreeWeight<ParsimonyScore> weight(dag);

  size_t score =
      weight.ComputeWeightBelow(dag.GetDAG().GetRoot(), {});

  assert_equal(score, expected_score, "Parsimony score");
}

static void test_subtree_weight(std::string_view path, size_t expected_score) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_subtree_weight(dag, expected_score);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_subtree_weight("data/testcase/full_dag.pb.gz", 75); },
              "Subtree weight: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_subtree_weight("data/testcase1/full_dag.pb.gz", 75); },
              "Subtree weight: testcase1"});
