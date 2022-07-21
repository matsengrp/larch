#include "subtree_weight.hpp"
#include "parsimony_score.hpp"

#include <iostream>

#include "test_common.hpp"

#include "synthetic_dags.hpp"

static void test_subtree_weight() {
  MADAG dag = MakeSyntheticDAG();

  dag.GetEdgeMutations() = dag.ComputeEdgeMutations(dag.GetReferenceSequence());

  SubtreeWeight<ParsimonyScore::Weight> weight{dag};

  size_t parsimony =
      weight.ComputeWeightBelow(dag.GetDAG().GetRoot(), ParsimonyScore{});

  assert_equal(parsimony, size_t(44), "Parsimony score");
}

[[maybe_unused]] static const auto test_added =
    add_test({[] { test_subtree_weight(); }, "Subtree weight"});
