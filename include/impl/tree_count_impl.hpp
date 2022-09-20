#include <algorithm>

TreeCount::Weight TreeCount::ComputeLeaf(const MADAG& dag, NodeId node_id) {
  return 1;
}

TreeCount::Weight TreeCount::ComputeEdge(const MADAG& dag, EdgeId edge_id) {
  /* This doesn't matter because AboveNode ignores edge weight */
  return 1;
}

std::pair<TreeCount::Weight, std::vector<size_t>> TreeCount::WithinCladeAccumOptimum(std::vector<TreeCount::Weight> inweights) {
    std::vector<size_t> indices;
    std::iota(indices.begin(), indices.end(), 0);
    return {std::accumulate(inweights.begin(), inweights.end(), 0, std::plus<TreeCount::Weight>()), indices};
}

TreeCount::Weight TreeCount::BetweenClades(std::vector<TreeCount::Weight> inweights) {
    return std::accumulate(inweights.begin(), inweights.end(), 1, std::multiplies<TreeCount::Weight>());
}

TreeCount::Weight TreeCount::AboveNode(TreeCount::Weight edgeweight, TreeCount::Weight childnodeweight) {
    return  childnodeweight;
}
