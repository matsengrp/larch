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
    dag.RecomputeEdgeMutations();
  }

  SubtreeWeight<WeightAccumulator<ParsimonyScore>> parsimonycount(dag);

  SubtreeWeight<ParsimonyScore> parsimonyscore(dag);

  Weight score = parsimonycount.ComputeWeightBelow(dag.GetDAG().GetRoot(), {});

  auto opt_pair = parsimonyscore.ComputeOptimalSubtreeCountBelow(dag.GetDAG().GetRoot(), {});

  assert_equal(score.GetWeights().begin()->first, opt_pair.first, "optimal weight doesn't match");
  assert_equal(score.GetWeights().begin()->second, opt_pair.second, "optimal weight count doesn't match");

  /* std::cout << "parsimony score counts\n"; */
  /* std::cout << "score   |   count\n"; */
  /* for (auto& scorepair : score.GetWeights()) { */
  /*   std::cout << scorepair.first << " | " << scorepair.second << "\n"; */
  /* } */
  assert_equal(score, expected_score, "unexpected parsimony counts");
}

static void test_weight_accum(std::string_view path, Weight expected_score) {
  MADAG dag = LoadDAGFromProtobuf(path);
  test_weight_accum(dag, expected_score);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_weight_accum("data/testcase/full_dag.pb.gz", Weight({{78, 211},
                                                                          {77, 206},
                                                                          {79, 143},
                                                                          {76, 106},
                                                                          {80, 79},
                                                                          {81, 27},
                                                                          {75, 23},
                                                                          {82, 11},
                                                                          {83, 9},
                                                                          {84, 3}},
                                                                         {}));
              },
              "Parsimony Counting: testcase"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] {
                test_weight_accum("data/testcase1/full_dag.pb.gz",
                                  Weight({{78, 3}, {79, 2}, {76, 1}, {75, 1}}, {}));
              },
              "Parsimony Counting: testcase1"});
