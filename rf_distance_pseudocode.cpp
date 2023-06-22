computing the summed RF distance over all trees in the DAG.
This routine requires 4 traversals of the nodes in DAG to calculate the rf distance sums: 
- one postorder traversal to assign subtree counts below each node, and
- one in a preorder traversal that uses the subtree counts and computes above-tree counts
- a loop over the nodes (in any order) accumulates the two types of counts and saves them in an accumulation based on unique leafsets
- finally, a postorder traversal to calculate the sums of RF distances between trees in a dag and all the trees in a reference, using the methods provided by a SubtreeWeight<SumRFDistance, DAGTYPE> object.

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

//Create a WeightOps for computing sum RF distances to the provided reference DAG:
template <typename REFDAG>
struct SumRFDistance {
  using Weight = ArbitraryInt;
  ArbitraryInt num_trees_in_dag;
  std::unordered_map<LeafSet, ArbitraryInt> leafset_to_full_treecount;
  ArbitraryInt shift_sum;

  SumRFDistance(REFDAG reference_dag) {
    SubtreeWeight<TreeCount, MADAG> below_tree_counts{reference_dag};
    std::vector<ArbitraryInt> above_tree_counts;
    above_tree_counts.resize(reference_dag.GetNodesCount());
    num_trees_in_dag = below_tree_counts(reference_dag.GetRoot());
    auto above_root_node = ComputeAboveTreeCount(reference_dag.GetRoot(), above_tree_counts, below_tree_counts);

    // create a list of unique (topologically) nodes in the DAG, and accumulate above_tree_counts[n]*below_tree_counts[n] by adding over all n with identical clade sets
    for (auto node: reference_dag.GetResult().GetNodes()) {
      if (not node.IsUA()) {
        leafset_to_full_treecount[reference_dag.GetResultNodeLabels().at(node.GetId().value).GetLeafSet()] += above_tree_counts[node.GetId().value]*below_tree_counts[node];
      }
    }

    // sum all of the values in leafset_to_full_treecount
    shift_sum = std::accumulate(leafset_to_full_treecount.begin(), leafset_to_full_treecount.end(),
    0, [](std::unordered_map<LeafSet, ArbitraryInt> v1, std::unordered_map<LeafSet, ArbitraryInt> v2){return v1.second + v2.second;}
    );

  }

  template <typename DAG>
  static Weight ComputeLeaf(DAG dag, NodeId node_id){
    return 0;
  }

  template <typename DAG>
  static Weight ComputeEdge(DAG dag, EdgeId edge_id){
    auto edge = dag.GetEdge(edge_id);
    //clade should be a LeafSet, the key type in the mapping
    //leafset_to_full_treecount:
    auto clade = edge.parent.GetClade(edge.clade);
    auto record = leafset_to_full_treecount.find(clade);
    if (record == leafset_to_full_treecount.end()) {
      return num_trees_in_dag;
    } else {
      return num_trees_in_dag - (2 * record->second);
    }
  }
  /*
   * Given a vector of weights for edges below a clade, compute the minimum
   * weight of them all, and return that minimum weight, and a vector
   * containing the indices of all elements of the passed vector that achieve
   * that minimum
   * TODO the following three methods are exactly the same as the implementation in
   * parsimony_score_impl.hpp, etc. Maybe we can avoid repeating everywhere?
   */
  inline static std::pair<Weight, std::vector<size_t>> WithinCladeAccumOptimum(
      const std::vector<Weight>&){
    Weight optimal_weight = std::numeric_limits<size_t>::max();
    std::vector<size_t> optimal_indices;
    size_t inweight_idx = 0;
    for (auto weight : inweights) {
      if (weight < optimal_weight) {
        optimal_weight = weight;
        optimal_indices.clear();
        optimal_indices.push_back(inweight_idx);
      } else if (weight == optimal_weight) {
        optimal_indices.push_back(inweight_idx);
      }
      inweight_idx++;
    }
    return {optimal_weight, optimal_indices};
  }
  /*
   * Given a vector of weights, one for each child clade, aggregate them
   */
  inline static Weight BetweenClades(const std::vector<Weight>&){
    return std::accumulate(inweights.begin(), inweights.end(),
                           static_cast<ParsimonyScore::Weight>(0));
  }
  inline static Weight AboveNode(Weight edgeweight, Weight childnodeweight){
    return edgeweight + childnodeweight;
  }
  //Will: I think the above WeightOps methods achieve all that's written in the following TODO,
  //although the resulting values from
  //using this WeightOps with SubtreeWeight will need to be have shift_sum added to them
  //to get correct summed RF distances.
  /* //TODO: create a SubtreeWeight that has: */
  /* leaf_func = lambda leaf_node: shift_sum */
  /* edge_func = lambda edge: num_trees_in_dag - 2*leaf_set_to_full_treecount[reference_dag.GetResultNodeLabels.at(edge.GetChild().GetId().value).GetLeafSet()] */
  /* clade_func = lambda clade: min(edge_func(e) for e in clade) */
  /* node_func = lambda internal_node: sum(clade_func(c) - shift_sum for c in internal_node.GetClades()) + shift_sum */
} 

//Create a WeightOps for computing RF distances to the provided reference tree:
template <typename REFDAG>
struct RFDistance : SumRFDistance{
  RFDistance(REFDAG reference_dag) {
    assert reference_dag.IsTree();
    //now behave exactly like SumRFDistance
  }
} 
