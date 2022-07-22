#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>
#include <string_view>

#include "test_common.hpp"

#include "synthetic_dags.hpp"
#include "dag_loader.hpp"

static void test_subtree_weight(MADAG& dag, size_t expected_score) {
  if (dag.GetEdgeMutations().empty()) {
    dag.GetEdgeMutations() = dag.ComputeEdgeMutations(dag.GetReferenceSequence());
  }

  SubtreeWeight<ParsimonyScore::Weight, ParsimonyScore> weight{dag};

  size_t score = weight.ComputeWeightBelow(dag.GetDAG().GetRoot(), ParsimonyScore{});

  assert_equal(score, expected_score, "Parsimony score");
}

static void test_subtree_weight(size_t expected_score) {
  MADAG dag = MakeSyntheticDAG();
  test_subtree_weight(dag, expected_score);
}

static void test_subtree_weight(std::string_view path, size_t expected_score) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_subtree_weight(dag, expected_score);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_subtree_weight(34); }, "Subtree weight: synthetic"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_subtree_weight("data/test_5_trees/tree_0.pb.gz", 174); },
              "Subtree weight: test_5_trees"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_subtree_weight("data/testcase2/full_dag.pb.gz", 174); },
              "Subtree weight: testcase2"});

[[maybe_unused]] static const auto test_added3 =
    add_test({[] { test_subtree_weight("data/testcaseref/tree_0_newref.pb.gz", 650); },
              "Subtree weight: testcaseref"});

[[maybe_unused]] static const auto test_added4 =
    add_test({[] { test_subtree_weight("data/20D_from_fasta/full_dag.pb.gz", 11149); },
              "Subtree weight: 20D_from_fasta"});
