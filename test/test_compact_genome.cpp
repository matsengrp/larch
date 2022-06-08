#include "compact_genome.hpp"

#include "test_common.hpp"
#include "dag_loader.hpp"
#include "merge.hpp"

static void test_compact_genome(std::string_view pb_path, std::string_view json_path) {
  MADAG json_dag = LoadDAGFromJson(json_path);
  MADAG pb_dag = LoadDAGFromProtobuf(pb_path);
  assert_equal(json_dag.reference_sequence, pb_dag.reference_sequence,
               "Reference sequence");

  // std::vector<CompactGenome> pb_cgs1 =
  //     pb_dag.ComputeCompactGenomes(pb_dag.reference_sequence);
  // std::vector<CompactGenome> pb_cgs2 =
  //     pb_dag.ComputeCompactGenomesDAG(pb_dag.reference_sequence);
  // assert_equal(pb_cgs1, pb_cgs2, "Compact genomes tree/DAG");

  // std::vector<CompactGenome> json_cgs = LoadCompactGenomesJson(json_path);
  // assert_equal(pb_cgs2, json_cgs, "Compact genomes");
}

static void test_edge_mutations(std::string_view path) {
  MADAG tree = LoadDAGFromProtobuf(path);
  Merge merge{tree.reference_sequence};
  std::vector<std::reference_wrapper<const MADAG>> tree_refs;
  tree_refs.push_back(tree);

  merge.AddTrees(tree_refs, false);

  std::vector<EdgeMutations> computed_mutations = merge.ComputeResultEdgeMutations();

  size_t failed = 0;
  for (Edge edge : merge.GetResult().GetEdges()) {
    const EdgeMutations& mutations = computed_mutations.at(edge.GetId().value);
    if (mutations != tree.edge_mutations.at(edge.GetId().value)) {
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

[[maybe_unused]] static const auto test_added4 =
    add_test({[] {
                test_compact_genome("data/test_5_trees/full_dag.pb.gz",
                                    "data/test_5_trees/full_dag.json.gz");
              },
              "Compact genomes: test_5_trees"});

[[maybe_unused]] static const auto test_added5 =
    add_test({[] {
                test_compact_genome("data/testcase2/full_dag.pb.gz",
                                    "data/testcase2/full_dag.json.gz");
              },
              "Compact genomes: testcase2"});

[[maybe_unused]] static const auto test_added6 =
    add_test({[] {
                test_compact_genome("data/20D_from_fasta/20D_full_dag.pb.gz",
                                    "data/20D_from_fasta/20D_full_dag.json.gz");
              },
              "Compact genomes: 20D_from_fasta"});
