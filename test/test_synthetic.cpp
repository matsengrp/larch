#include "synthetic_dags.hpp"

#include <iostream>

#include "test_common.hpp"

static void test_compact_genome() {
  MADAG dag = MakeSyntheticDAG();

  std::vector<EdgeMutations> computed_mutations =
      dag.ComputeEdgeMutations(dag.reference_sequence);

  Assert(dag.compact_genomes.size() == dag.dag.GetNodesCount());
  Assert(computed_mutations.size() == dag.dag.GetEdgesCount());

  std::vector<CompactGenome> computed_cgs =
      dag.ComputeCompactGenomes(dag.reference_sequence);

  assert_equal(computed_cgs, dag.compact_genomes, "Compact genomes");
}

[[maybe_unused]] static const auto test_added =
    add_test({[] { test_compact_genome(); }, "Synthetic: compact genomes"});
