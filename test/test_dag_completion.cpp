#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/dag_loader.hpp"

using Node = MutableMADAG::NodeView;
using Edge = MutableMADAG::EdgeView;

[[maybe_unused]] void dag_info(MADAGStorage& dag_storage) {
  auto dag = dag_storage.View();

  std::cout << std::endl << "=== DAG_INFO [begin] ===" << std::endl;
  std::cout << "=== NODES: " << dag.GetNodesCount() << " ===" << std::endl;
  for (auto node : dag.GetNodes()) {
    std::cout << "MADAG::Node: " << node.GetId() << " "
              << (node.IsLeaf() ? "LEAF" : "NON-LEAF") << std::endl;
    std::cout << "MADAG::CompactGenome: " << node.GetCompactGenome().ToString()
              << std::endl;
    std::cout << "MADAG::LeafSet: " << node.GetLeafsBelow() << std::endl;
    std::cout << "MADAG::Children: " << node.GetChildren() << std::endl;
    std::cout << "MADAG::Clades: " << node.GetClades() << std::endl;
  }
  std::cout << "=== EDGES: " << dag.GetEdgesCount() << " ===" << std::endl;
  for (auto edge : dag.GetEdges()) {
    std::cout << "MADAG::Edge: " << edge.GetId() << " [" << edge.GetParent() << " -> "
              << edge.GetChild() << " | " << edge.GetClade() << "]" << std::endl;
  }
  std::cout << "=== DAG_INFO [end] ===" << std::endl << std::endl;
}

[[maybe_unused]] std::set<std::tuple<NodeId, NodeId, CladeIdx>> dag_make_node_pair_map(
    MADAGStorage& dag_storage) {
  auto dag = dag_storage.View();
  std::set<std::tuple<NodeId, NodeId, CladeIdx>> node_pairs;
  for (auto edge : dag.GetEdges()) {
    node_pairs.insert({edge.GetParent(), edge.GetChild(), edge.GetClade()});
  }
  return node_pairs;
}

[[maybe_unused]] bool dag_compare_topologies(MADAGStorage& lhs_storage,
                                             MADAGStorage& rhs_storage) {
  auto lhs_node_map = dag_make_node_pair_map(lhs_storage);
  auto rhs_node_map = dag_make_node_pair_map(rhs_storage);
  return lhs_node_map == rhs_node_map;
}

[[maybe_unused]] void dag_make_complete_rootsplits(MADAGStorage& dag_storage) {
  auto dag = dag_storage.View();
  auto clade_union_map = dag.BuildCladeUnionMap();
  size_t leaf_count = 0;
  std::set<NodeId>* rootsplits = nullptr;
  for (auto& [clade_union, node_ids] : clade_union_map) {
    if (clade_union.size() > leaf_count) {
      rootsplits = &node_ids;
    }
    leaf_count = std::max(leaf_count, clade_union.size());
  }
  for (auto node_id : *rootsplits) {
    dag.AppendEdge(dag.GetRoot().GetId(), node_id, {0});
  }
}

[[maybe_unused]] void test_dag_completion() {
  auto dag_storage = MakeSampleDAG();
  auto dag = dag_storage.View();
  auto dag_storage_2 = MakeSampleDAG();

  dag.GetRoot().CalculateLeafsBelow();
  // dag_info(dag_storage);
  assert_true(dag_compare_topologies(dag_storage, dag_storage_2),
              "DAGs are not equal before removing edges.");

  dag.ClearConnections();
  // dag_info(dag_storage);
  assert_false(dag_compare_topologies(dag_storage, dag_storage_2),
               "DAGs are incorrectly equal after removing edges.");

  dag_make_complete_rootsplits(dag_storage);
  dag.MakeComplete();
  // dag_info(dag_storage);
  assert_true(dag_compare_topologies(dag_storage, dag_storage_2),
              "DAGs are not equal after recompleting the DAG.");
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_dag_completion(); }, "DAG Completion"});
