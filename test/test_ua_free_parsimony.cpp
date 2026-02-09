#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/ua_free_parsimony_score.hpp"
#include "larch/subtree/parsimony_score.hpp"

#include <iostream>
#include <string_view>
#include <string>

#include "test_common.hpp"

#include "larch/dag_loader.hpp"

static void test_ua_free_parsimony_score(MADAG dag) {
  SubtreeWeight<UAFreeParsimonyScore, MADAG> weight_ua{dag};
  SubtreeWeight<ParsimonyScore, MADAG> weight{dag};

  UAFreeParsimonyScore::Weight score_ua =
      weight_ua.ComputeWeightBelow(dag.GetRoot(), {});
  ParsimonyScore::Weight score = weight.ComputeWeightBelow(dag.GetRoot(), {});
  for (auto child : dag.GetRoot().GetChildren()) {
    score -= child.GetEdgeMutations().size();
  }
  TestAssert(score_ua == score);
}

static void test_ua_free_parsimony_score(std::string_view path,
                                         std::string_view refseq_path) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage dag = LoadTreeFromProtobuf(path, reference_sequence);
  dag.View().RecomputeCompactGenomes(true);
  test_ua_free_parsimony_score(dag.View());
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_ua_free_parsimony_score("data/seedtree/seedtree.pb.gz",
                                             "data/seedtree/refseq.txt.gz");
              },
              "UA Free parsimony: testcase"});
