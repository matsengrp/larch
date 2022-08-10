#include "compact_genome.hpp"

#include "test_common.hpp"
#include "dag_loader.hpp"
#include "merge.hpp"

[[maybe_unused]] static void test_edge_mutations(std::string_view path) {
  MADAG dag = LoadDAGFromProtobuf(path);
  dag.GetCompactGenomes() = dag.ComputeCompactGenomes(dag.GetReferenceSequence());

  std::vector<EdgeMutations> computed_mutations =
      dag.ComputeEdgeMutations(dag.GetReferenceSequence());

  Assert(dag.GetCompactGenomes().size() == dag.GetDAG().GetNodesCount());
  Assert(computed_mutations.size() == dag.GetDAG().GetEdgesCount());

  std::vector<CompactGenome> computed_cgs =
      dag.ComputeCompactGenomes(dag.GetReferenceSequence());

  assert_equal(computed_cgs, dag.GetCompactGenomes(), "Compact genomes");

  assert_equal(dag.ComputeEdgeMutations(dag.GetReferenceSequence()),
               dag.GetEdgeMutations(), "Edge mutations");
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_edge_mutations("data/test_5_trees/tree_0.pb.gz"); },
              "Compact genomes: test_5_trees"});
