#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>
#include <string_view>

#include "test_common.hpp"

#include "dag_loader.hpp"

static void test_sample_tree(MADAG& dag) {
  if (dag.GetEdgeMutations().empty()) {
    dag.RecomputeEdgeMutations();
  }

  SubtreeWeight<ParsimonyScore> weight(dag);

  MADAG result = weight.SampleTree({}).first;

  assert_true(result.GetDAG().IsTree(), "Tree");
}

static void test_sample_tree(std::string_view path) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_sample_tree(dag);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_sample_tree("data/testcase/full_dag.pb.gz"); },
              "Sample tree: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_sample_tree("data/testcase1/full_dag.pb.gz"); },
              "Sample tree: testcase1"});
