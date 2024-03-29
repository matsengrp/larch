#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"

#include <iostream>
#include <string_view>

#include "test_common.hpp"

#include "larch/dag_loader.hpp"

static void test_subtree_weight(MADAG dag, size_t expected_score) {
  SubtreeWeight<BinaryParsimonyScore, MADAG> weight(dag);

  size_t score = weight.ComputeWeightBelow(dag.GetRoot(), {});

  TestAssert(score == expected_score);
}

static void test_subtree_weight(std::string_view path, size_t expected_score) {
  MADAGStorage dag = LoadDAGFromProtobuf(path);
  auto dag_view = dag.View();
  dag_view.RecomputeCompactGenomes(true);
  dag_view.SampleIdsFromCG(true);
  test_subtree_weight(dag_view, expected_score);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_subtree_weight("data/testcase/full_dag.pb.gz", 75); },
              "Subtree weight: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_subtree_weight("data/testcase1/full_dag.pb.gz", 75); },
              "Subtree weight: testcase1"});
