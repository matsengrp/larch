#include "compact_genome.hpp"

#include "test_common.hpp"
#include "dag_loader.hpp"
#include "merge.hpp"

static void do_test(std::string_view path) {
  std::string reference_sequence;
  std::vector<std::vector<Mutations>> edge_mutations;
  std::vector<DAG> trees;
  edge_mutations.push_back({});
  trees.push_back(LoadDAGFromProtobuf(path, reference_sequence, edge_mutations.at(0)));
  Merge merge{reference_sequence};
  std::vector<std::reference_wrapper<const DAG>> tree_refs{trees.begin(), trees.end()};
  merge.AddTrees(tree_refs, edge_mutations, false);

  std::vector<Mutations> computed_mutations = merge.ComputeResultEdgeMutations();

  size_t failed = 0;
  for (Edge edge : merge.GetResult().GetEdges()) {
    const Mutations& mutations = computed_mutations.at(edge.GetId().value);
    if (mutations != edge_mutations.at(0).at(edge.GetId().value)) {
      ++failed;
    }
  }
  if (failed > 0) {
    std::cerr << "Failed: " << failed << "\n";
  }
  Assert(failed == 0);
}

static void test_5_trees() { do_test("data/test_5_trees/tree_0.pb.gz"); }

static void test_800_trees() { do_test("data/20D_from_fasta/20D_full_dag.pb.gz"); }

[[maybe_unused]] static const auto test_added0 =
    add_test({test_5_trees, "Compact genome: 5 trees"});

[[maybe_unused]] static const auto test_added2 =
    add_test({test_800_trees, "Compact genome: 800 trees"});
