#include "compact_genome.hpp"

#include "test_common.hpp"
#include "dag_loader.hpp"
#include "merge.hpp"

[[maybe_unused]] static void test_edge_mutations(std::string_view path) {
  MADAG dag = LoadDAGFromProtobuf(path);
  dag.compact_genomes = dag.ComputeCompactGenomes(dag.reference_sequence);

  std::vector<EdgeMutations> computed_mutations =
      dag.ComputeEdgeMutations(dag.reference_sequence);

  Assert(dag.compact_genomes.size() == dag.dag.GetNodesCount());
  Assert(computed_mutations.size() == dag.dag.GetEdgesCount());

  std::vector<CompactGenome> computed_cgs =
      dag.ComputeCompactGenomes(dag.reference_sequence);

  assert_equal(computed_cgs, dag.compact_genomes, "Compact genomes");

  assert_equal(dag.ComputeEdgeMutations(dag.reference_sequence), dag.edge_mutations,
               "Edge mutations");
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_edge_mutations("data/test_5_trees/tree_0.pb.gz"); },
              "Compact genomes: test_5_trees"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_edge_mutations("data/testcase2/full_dag.pb.gz"); },
              "Compact genomes: testcase2"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_edge_mutations("data/testcaseref/tree_0_newref.pb.gz"); },
              "Compact genomes: testcaseref"});

[[maybe_unused]] static const auto test_added3 =
    add_test({[] { test_edge_mutations("data/20D_from_fasta/full_dag.pb.gz"); },
              "Compact genomes: 20D_from_fasta"});
