#include "test_common.hpp"
#include "larch/dag_loader.hpp"

// set<CompactGenome> is just to express the spirit of what we want.
using CladeUnion = std::set<CompactGenome>;

CladeUnion Node::GetCladeUnion() {
  // return this node's LeafSet, but flattened, so that it's just
  // a set/ordered vector of CompactGenomes
  CladeUnion clade_union;

  return clade_union;
}

template <typename DAGType, typename NodeType>
[[maybe_unused]] CompleteDAG(DAGType dag) {
  std::map<CladeUnion, std::vector<NodeType>> target_nodes_map;
  for (auto node : dag.GetNodes()) {
    auto node_vec = target_nodes_map[node.GetCladeUnion()];
    node_vec.insert(node);
  }

  for (auto node : dag.GetNodes()) {
    // assuming that LeafSet indices align with clade indices
    LeafSet clades = node.GetLeafSet();
    for (size_t clade_idx = 0; i < clades.size(); i++) {
      if (auto possible_children = target_nodes_map.find(
              clades[i]; possible_children != target_nodes_map.end())) {
        for (auto child_node : possible_children->second) {
          // of course, we don't want to duplicate edges that are already in
          // the DAG, but I assume that AddEdge handles this?
          dag.AddEdge(/*parent*/ node, /*child*/ child_node, /*clade index*/ clade_idx);
        }
      }
    }
  }
}

[[maybe_unused]] void test_dag_completion() {
  std::cout << "test_dag_completion [begin]" << std::endl;

  std::cout << "test_dag_completion [end]" << std::endl;
}

[[maybe_unused]] static const auto test_added0 = add_test(
    {[] { test_edge_mutations("data/test_5_trees/tree_0.pb.gz"); }, "DAG Completion"});
