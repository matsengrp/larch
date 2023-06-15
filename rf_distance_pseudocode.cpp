computing the summed RF distance over all trees in the DAG.
This routine requires 4 traversals of the nodes in DAG to calculate the rf distance sums: 
- one postorder traversal to assign subtree counts below each node, and
- one in a preorder traversal that uses the subtree counts and computes above-tree counts
- a loop over the nodes (in any order) accumulates the two types of counts and saves them in an accumulation based on unique leafsets
- finally, a postorder traversal to calculate the 

template <typename subtreeWeight>
ArbitraryInt ComputeAboveTreeCount(Node node, std::vector<ArbitraryInt>& above_tree_counts, subtreeWeight cached_below_tree_counts) {
  /* This routine takes an argument that is a cached subtreeWeight `cached_below_tree_counts'. If this subtreeWeight is the TreeCount weight (which returns the number of subtrees below each node), then this routine computes the number of trees "above" a given node. A tree is "above" the node if it contains the node, and taking the graph union with a subtree below that node yields a tree on the full leaf set that belongs in the DAG.
  */

  if (above_tree_counts.at(node.GetId().value) > 0) {
    return above_tree_counts.at(node.GetId().value);
  } else if (node.IsUANode()) {
    above_tree_counts.at(node.GetId().value) = 1;
    return 1;
  }

  ArbitraryInt above = 0;
  for (auto parent_edge: node.GetParents()) {
    auto parent_node = parent_edge.GetParent();
    auto current_clade = parent_edge.GetClade();
    ArbitraryInt below_parent = 1;

    for (auto clade: parent_node.GetClades()) {
      if (clade != current_clade) {
        // sum up the list of values {cached_below_tree_counts[edge] : edge in clade}
        auto below_clade = ranges::views::accumulate(parent_node.GetClade(clade) | ranges::views::transform::([](Edge &edge) {return cached_below_tree_counts(edge.GetChild().GetId().value);}), 0);
        // multiply these sums across all the (other) clades below the current parent node
        below_parent *= below_clade;
      }
    }
    // sum all values over parents for current node
    above += ComputeAboveTreeCount(parent_node, above_tree_counts, cached_below_tree_counts) * below_parent;
  }
  above_tree_counts.at(node.getId().value) = above;
  return above;
}

template <typename MERGE>
ArbitraryInt ComputeSumRFDistance(MERGE merge, DAG dag) {

  SubtreeWeight<TreeCount, MADAG> below_tree_counts{merge.GetResult()};
  std::vector<ArbitraryInt> above_tree_counts;
  above_tree_counts.resize(merge.GetResult().GetNodesCount());
  auto num_trees_in_dag = below_tree_counts(merge.GetResult().GetRoot());
  auto above_root_node = ComputeAboveTreeCount(merge.GetResult().GetRoot(), above_tree_counts, below_tree_counts);

  // create a list of unique (topologically) nodes in the DAG, and accumulate above_tree_counts[n]*below_tree_counts[n] by adding over all n with identical clade sets
  std::unordered_map<LeafSet, ArbitraryInt> leafset_to_full_treecount;
  for (auto node: merge.GetResult().GetNodes()) {
    if (not node.IsUA()) {
      leafset_to_full_treecount[merge.GetResultNodeLabels().at(node.GetId().value).GetLeafSet()] += above_tree_counts[node.GetId().value]*below_tree_counts[node];
    }
  }

  // sum all of the values in leafset_to_full_treecount
  ArbitraryInt shift_sum = std::accumulate(leafset_to_full_treecount.begin(), leafset_to_full_treecount.end(),
  0, [](std::unordered_map<LeafSet, ArbitraryInt> v1, std::unordered_map<LeafSet, ArbitraryInt> v2){return v1.second + v2.second;}
  );

  //TODO: create a SubtreeWeight that has:
  leaf_func = lambda leaf_node: shift_sum
  edge_func = lambda edge: num_trees_in_dag - 2*leaf_set_to_full_treecount[merge.GetResultNodeLabels.at(edge.GetChild().GetId().value).GetLeafSet()]
  clade_func = lambda clade: min(edge_func(e) for e in clade)
  node_func = lambda internal_node: sum(clade_func(c) - shift_sum for c in internal_node.GetClades()) + shift_sum
}


