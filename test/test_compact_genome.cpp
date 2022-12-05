#include "larch/madag/compact_genome.hpp"

#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/merge/merge.hpp"

[[maybe_unused]] static void test_edge_mutations(std::string_view path) {
  MADAGStorage dag_storage = LoadDAGFromProtobuf(path);
  MutableMADAG dag = dag_storage.View();
  using Node = MutableMADAG::NodeView;
  using Edge = MutableMADAG::EdgeView;

  std::vector<EdgeMutations> loaded_edge_mutatons;
  for (Edge edge : dag.GetEdges()) {
    loaded_edge_mutatons.emplace_back(edge.GetEdgeMutations().Copy());
  }

  dag.RecomputeCompactGenomes();
  std::vector<CompactGenome> computed_cgs;
  for (Node node : dag.GetNodes()) {
    computed_cgs.emplace_back(node.GetCompactGenome().Copy());
  }
  for (Edge edge : dag.GetEdges()) {
    edge.SetEdgeMutations({});
  }
  dag.RecomputeEdgeMutations();

  size_t index = 0;
  for (Edge edge : dag.GetEdges()) {
    assert_equal(edge.GetEdgeMutations(), loaded_edge_mutatons.at(index++),
                 "Edge mutations");
  }

  for (Node node : dag.GetNodes()) {
    node.SetCompactGenome({});
  }
  dag.RecomputeCompactGenomes();
  index = 0;
  for (Node node : dag.GetNodes()) {
    assert_equal(node.GetCompactGenome(), computed_cgs.at(index++), "Compact genomes");
  }
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_edge_mutations("data/test_5_trees/tree_0.pb.gz"); },
              "Compact genomes: test_5_trees"});
