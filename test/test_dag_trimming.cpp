#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>
#include <string_view>

#include "test_common.hpp"

#include "dag_loader.hpp"

static void test_dag_trimming(MADAG& dag, size_t expected_score) {
  if (dag.GetEdgeMutations().empty()) {
    dag.GetEdgeMutations() = dag.ComputeEdgeMutations(dag.GetReferenceSequence());
  }

  SubtreeWeight<ParsimonyScore> weight(dag.GetDAG());

  size_t score = weight.ComputeWeightBelow(dag.GetDAG().GetRoot(), ParsimonyScore{dag});

  assert_equal(score, expected_score, "Parsimony score");
}

static void test_dag_trimming(std::string_view path, size_t expected_score) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_dag_trimming(dag, expected_score);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_dag_trimming("data/testcase/full_dag.pb.gz", 75); },
              "DAG trimming: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_dag_trimming("data/testcase1/full_dag.pb.gz", 75); },
              "DAG trimming: testcase1"});
