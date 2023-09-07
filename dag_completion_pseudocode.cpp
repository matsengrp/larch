/* Describes "history DAG completion", which means that all allowed edges are
 * added between the nodes already in the hDAG. An edge is allowed if the clade
 * it descends from matches the clade union of the target node. So, we start
 * with a mapping from clade unions (i.e. flattened LeafSets) to vectors of
 * existing hDAG nodes, then iterate through nodes v in the hDAG and for each
 * child clade c of the node, and for each other node v' in the dag whose clade
 * union matches that child clade, add an edge from v to v' descending from the
 * child clade c.
 *
 *
 * In practice there aren't usually very many edges to add, and completing the
 * hDAG often finds more MP trees, or sometimes even trees with better
 * parsimony scores. So, we could occasionally complete, then trim, the hDAG
 * during larch-usher optimization.
 *
 * I'm not sure where this should be implemented (perhaps as a method on the
 * merge object?) so I'll just write the pseudocode as a standalone function that
 * modifies the provided hDAG
 */

// set<CompactGenome> is just to express the spirit of what we want.
using CladeUnion = std::set<CompactGenome>;

template <typename NodeType>
CladeUnion GetCladeUnion(NodeType node) {
  // return this node's LeafSet, but flattened, so that it's just
  // a set/ordered vector of CompactGenomes
  CladeUnion clade_union;
  for (const auto compact_genome : node.GetCompactGenomes()) {
    clade_union.insert(compact_genomes);
  }
  return clade_union;
}

template <typename DAGType, NodeType>
void complete_dag(DAGTYPE dag) {
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
