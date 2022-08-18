#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>
#include <string_view>

#include "test_common.hpp"

#include "dag_loader.hpp"

static void test_sample_tree(MADAG& dag, size_t expected_score) {
  if (dag.GetEdgeMutations().empty()) {
    dag.GetEdgeMutations() = dag.ComputeEdgeMutations(dag.GetReferenceSequence());
  }

  SubtreeWeight<ParsimonyScore> weight(dag);

  MADAG result = weight.SampleTree({});

  assert_true(result.GetDAG().GetEdgesCount() > 0, "Edges");
}

static void test_sample_tree(std::string_view path, size_t expected_score) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_sample_tree(dag, expected_score);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_sample_tree("data/testcase/full_dag.pb.gz", 75); },
              "Sample tree: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_sample_tree("data/testcase1/full_dag.pb.gz", 75); },
              "Sample tree: testcase1"});
