#include "compact_genome.hpp"

#include "test_common.hpp"
#include "dag_loader.hpp"
#include "merge.hpp"

[[maybe_unused]] static void test_edge_mutations(std::string_view path) {
  MADAG tree = LoadDAGFromProtobuf(path);
  tree.compact_genomes = tree.ComputeCompactGenomesDAG(tree.reference_sequence);

  std::vector<EdgeMutations> computed_mutations =
      tree.ComputeEdgeMutations(tree.reference_sequence);

  Assert(tree.compact_genomes.size() == tree.dag.GetNodes().size());
  Assert(computed_mutations.size() == tree.dag.GetEdges().size());

  size_t failed = 0;
  for (EdgeId edge : tree.dag.GetEdges()) {
    const EdgeMutations& mutations = computed_mutations.at(edge.value);
    if (mutations != tree.edge_mutations.at(edge.value)) {
      ++failed;
    }
  }
  if (failed > 0) {
    std::cerr << "Failed: " << failed << "\n";
  }
  assert_equal(failed, size_t{0}, "Edge mutations");
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_edge_mutations("data/test_5_trees/tree_0.pb.gz"); },
              "Edge mutations: test_5_trees"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_edge_mutations("data/testcase2/full_dag.pb.gz"); },
              "Edge mutations: testcase2"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_edge_mutations("data/testcaseref/tree_1_newref.pb.gz"); },
              "Edge mutations: testcaseref"});

[[maybe_unused]] static const auto test_added3 =
    add_test({[] { test_edge_mutations("data/20D_from_fasta/20D_full_dag.pb.gz"); },
              "Edge mutations: 20D_from_fasta"});
