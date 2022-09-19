#include <algorithm>

TreeCount::Weight TreeCount::ComputeLeaf(const MADAG& dag, NodeId node_id) {
  return 1;
}

TreeCount::Weight TreeCount::ComputeEdge(const MADAG& dag, EdgeId edge_id) {
  /* This doesn't matter because AboveNode ignores edge weight */
  return 1;
}

std::pair<TreeCount::Weight, std::vector<size_t>> WithinCladeAccumOptimum(std::vector<TreeCount::Weight> inweights) {
    TreeCount::Weight result = 1;
    std::vector<size_t> optimal_indices;
    size_t inweight_idx = 0;
    for (auto weight : inweights) {
        result += weight;
        optimal_indices.push_back(inweight_idx);
        inweight_idx++;
    }
    return {result, optimal_indices};
}

TreeCount::Weight BetweenClades(std::vector<TreeCount::Weight> inweights) {
    return std::accumulate(inweights.begin(), inweights.end(), 1, std::multiplies<TreeCount::Weight>());
}

TreeCount::Weight AboveNode(TreeCount::Weight edgeweight, TreeCount::Weight childnodeweight) {
    return  childnodeweight;
}
