#include "subtree_weight.hpp"
#include "parsimony_score.hpp"
#include "weight_accumulator.hpp"

#include <iostream>
#include <string_view>
#include <string>

#include "test_common.hpp"

#include "dag_loader.hpp"

using Weight = typename WeightAccumulator<ParsimonyScore>::Weight;

static void test_weight_accum(MADAG& dag, Weight expected_score) {
  if (dag.GetEdgeMutations().empty()) {
    dag.GetEdgeMutations() = dag.ComputeEdgeMutations(dag.GetReferenceSequence());
  }

  SubtreeWeight<WeightAccumulator<ParsimonyScore>> parsimonycount(dag);

  Weight score = parsimonycount.ComputeWeightBelow(dag.GetDAG().GetRoot(), {});

  std::cout << "parsimony score counts";
  std::cout << "score   |   count";
  for (auto& scorepair : score.GetWeights()) {
      std::cout << scorepair.first << " | " << scorepair.second;
  }
  /* assert_equal(score, expected_score, "expected parsimony counts"); */
}

static void test_weight_accum(std::string_view path, Weight expected_score) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_weight_accum(dag, expected_score);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_weight_accum("data/testcase/full_dag.pb.gz", Weight({})); },
              "Subtree weight: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_weight_accum("data/testcase1/full_dag.pb.gz", Weight({})); },
              "Subtree weight: testcase1"});