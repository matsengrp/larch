#include <algorithm>

TreeCount::Weight TreeCount::ComputeLeaf(MADAG, NodeId) { return 1; }

TreeCount::Weight TreeCount::ComputeEdge(MADAG, EdgeId) {
  /* This doesn't matter because AboveNode ignores edge weight */
  return 1;
}

std::pair<TreeCount::Weight, std::vector<size_t>> TreeCount::WithinCladeAccumOptimum(
    std::vector<TreeCount::Weight> inweights) {
  std::vector<size_t> indices;
  std::iota(indices.begin(), indices.end(), 0);
  return {std::accumulate(inweights.begin(), inweights.end(),
                          static_cast<TreeCount::Weight>(0), std::plus<>()),
          indices};
}

TreeCount::Weight TreeCount::BetweenClades(
    const std::vector<TreeCount::Weight>& inweights) {
  return std::accumulate(inweights.begin(), inweights.end(),
                         static_cast<TreeCount::Weight>(1), std::multiplies<>());
}

TreeCount::Weight TreeCount::AboveNode(const TreeCount::Weight&,
                                       TreeCount::Weight childnodeweight) {
  return childnodeweight;
}
