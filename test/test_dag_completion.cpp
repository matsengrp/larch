#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/dag_loader.hpp"

using Node = MutableMADAG::NodeView;
using Edge = MutableMADAG::EdgeView;
// set<CompactGenome> is just to express the spirit of what we want.
using CladeUnion = std::set<CompactGenome>;
using CladeUnionMap = std::map<CladeUnion, std::vector<NodeId>>;

template <typename DAGType>
CladeUnionMap BuildCladeUnionMap(DAGType &dag) {
  CladeUnionMap clade_union_map;
  for (auto node : dag.GetNodes()) {
    std::cout << "Node: " << node.GetId() << std::endl;
    std::cout << "CompactGenome: " << node.GetCompactGenome().ToString() << std::endl;
    // auto node_vec = target_nodes_map[GetCladeUnion(node)];
    // node_vec.insert(node);
  }
  return clade_union_map;
}

// std::vector<std::vector<const SampleId*>> clades_union(
//     const std::vector<std::vector<const SampleId*>>& lhs,
//     const std::vector<std::vector<const SampleId*>>& rhs) {
//   std::vector<std::vector<const SampleId*>> result;

//   for (auto [lhs_clade, rhs_clade] : ranges::views::zip(lhs, rhs)) {
//     std::vector<const SampleId*> clade{lhs_clade};
//     clade.insert(clade.end(), rhs_clade.begin(), rhs_clade.end());
//     ranges::sort(clade);
//     ranges::unique(clade);
//     result.push_back(std::move(clade));
//   }

//   ranges::sort(result);
//   return result;
// }

CladeUnion GetCladeUnion(Node node) {
  std::cout << "GetCladeUnion: " << node.GetId() << std::endl;
  // return this node's LeafSet, but flattened, so that it's just
  // a set/ordered vector of CompactGenomes
  CladeUnion clade_union;
  std::cout << "CompactGenome: " << node.GetCompactGenome().ToString() << std::endl;
  // for (const auto &compact_genome : node.GetCompactGenome()) {
  //   clade_union.insert(compact_genome);
  //   std::cout << compact_genome.
  // }
  return clade_union;
}

template <typename DAGType>
[[maybe_unused]] void CompleteDAG(DAGType &dag) {
  auto clade_union_map = BuildCladeUnionMap(dag);

  for (auto node : dag.GetNodes()) {
    node.CalculateLeafsBelow();
  }

  for (auto node : dag.GetNodes()) {
    std::cout << "Node::FindChildren: " << node.GetId() << std::endl;
    // assuming that LeafSet indices align with clade indices
    auto clades = node.GetLeafsBelow();
    std::cout << "Node::Leafset: " << clades << std::endl;
    // for (size_t clade_idx = 0; clade_idx < clades.size(); clade_idx++) {
    //   if (auto possible_children = target_nodes_map.find(clades[clade_idx]);
    //       possible_children != target_nodes_map.end()) {
    //     for (auto child_node : possible_children->second) {
    //       // of course, we don't want to duplicate edges that are already in
    //       // the DAG, but I assume that AddEdge handles this?
    //       dag.AddEdge(/*parent*/ node, /*child*/ child_node,
    //                   /*clade index*/ clade_idx);
    //     }
    //   }
    // }
  }
}

[[maybe_unused]] void test_dag_completion() {
  std::cout << "test_dag_completion [begin]" << std::endl;

  auto dag = MakeSampleDAG();
  dag.CalculateLeafsBelow();
  CompleteDAG(dag);

  std::cout << "test_dag_completion [end]" << std::endl;
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_dag_completion(); }, "DAG Completion"});
