#include "synthetic_dags.hpp"

#include <iostream>

#include "test_common.hpp"

static void test_compact_genome() {
  MADAG dag = MakeSyntheticDAG();

  std::vector<EdgeMutations> computed_mutations = dag.ComputeEdgeMutations();

  Assert(dag.GetCompactGenomes().size() == dag.GetDAG().GetNodesCount());
  Assert(computed_mutations.size() == dag.GetDAG().GetEdgesCount());

  std::vector<CompactGenome> computed_cgs = dag.ComputeCompactGenomes();

  assert_equal(computed_cgs, dag.GetCompactGenomes(), "Compact genomes");
}

[[maybe_unused]] static const auto test_added =
    add_test({[] { test_compact_genome(); }, "Synthetic: compact genomes"});
