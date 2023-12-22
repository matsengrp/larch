#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/subtree/weight_accumulator.hpp"

#include <iostream>
#include <string_view>
#include <string>

#include "test_common.hpp"

#include "larch/dag_loader.hpp"

using Weight = typename WeightAccumulator<ParsimonyScore>::Weight;

static void test_weight_accum(MADAG dag, Weight expected_score) {
  SubtreeWeight<WeightAccumulator<ParsimonyScore>, MADAG> parsimonycount(dag);

  Weight score = parsimonycount.ComputeWeightBelow(dag.GetRoot(), {});

  /* std::cout << "parsimony score counts\n"; */
  /* std::cout << "score   |   count\n"; */
  /* for (auto& scorepair : score.GetWeights()) { */
  /*   std::cout << scorepair.first << " | " << scorepair.second << "\n"; */
  /* } */
  TestAssert(score == expected_score);
}

static void test_weight_accum(std::string_view path, Weight expected_score) {
  MADAGStorage dag = LoadDAGFromProtobuf(path);
  auto dag_view = dag.View();
  dag_view.RecomputeCompactGenomes(true);
  dag_view.SampleIdsFromCG(true);
  test_weight_accum(dag_view, std::move(expected_score));
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
