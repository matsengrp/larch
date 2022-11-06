#include <algorithm>

ParsimonyScore::Weight ParsimonyScore::ComputeLeaf(MADAG, NodeId) { return 0; }

ParsimonyScore::Weight ParsimonyScore::ComputeEdge(MADAG dag, EdgeId edge_id) {
  return dag.Get(edge_id).GetEdgeMutations().size();
}

std::pair<ParsimonyScore::Weight, std::vector<size_t>>
ParsimonyScore::WithinCladeAccumOptimum(
    const std::vector<ParsimonyScore::Weight>& inweights) {
  ParsimonyScore::Weight optimal_weight = std::numeric_limits<size_t>::max();
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

ParsimonyScore::Weight ParsimonyScore::BetweenClades(
    const std::vector<ParsimonyScore::Weight>& inweights) {
  return std::accumulate(inweights.begin(), inweights.end(),
                         static_cast<ParsimonyScore::Weight>(0));
}

ParsimonyScore::Weight ParsimonyScore::AboveNode(
    ParsimonyScore::Weight edgeweight, ParsimonyScore::Weight childnodeweight) {
  return edgeweight + childnodeweight;
}
