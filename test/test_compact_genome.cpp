#include "compact_genome.hpp"

#include "test_common.hpp"
#include "dag_loader.hpp"
#include "merge.hpp"

void test_compact_genome() {
  std::string reference_sequence;
  std::vector<std::vector<Mutations>> edge_mutations;
  std::vector<DAG> trees;
  edge_mutations.push_back({});
  trees.push_back(LoadDAGFromProtobuf("data/20D_from_fasta/20D_full_dag.pb.gz",
                                      reference_sequence, edge_mutations.at(0)));
  Merge merge{reference_sequence};
  std::vector<std::reference_wrapper<const DAG>> tree_refs{trees.begin(), trees.end()};
  merge.AddTrees(tree_refs, edge_mutations, false);

  std::vector<Mutations> computed_mutations = merge.ComputeResultEdgeMutations();

  for (Edge edge : merge.GetResult().GetEdges()) {
    const Mutations& mutations = computed_mutations.at(edge.GetId().value);
    Assert(mutations == edge_mutations.at(0).at(edge.GetId().value));
  }
}

[[maybe_unused]] static const auto test_added =
    add_test({test_compact_genome, "Compact genome"});
