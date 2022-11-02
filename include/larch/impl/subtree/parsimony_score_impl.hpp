#include <algorithm>

ParsimonyScore::Weight ParsimonyScore::ComputeLeaf(const MADAG&, NodeId) { return 0; }

ParsimonyScore::Weight ParsimonyScore::ComputeEdge(const MADAG& dag, EdgeId edge_id) {
  return dag.GetEdgeMutations(edge_id).size();
}

std::pair<ParsimonyScore::Weight, std::vector<size_t>>
ParsimonyScore::WithinCladeAccumOptimum(std::vector<ParsimonyScore::Weight> inweights) {
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
    const std::vector<ParsimonyScore::Weight>& inweights) const {
  return std::accumulate(inweights.begin(), inweights.end(),
                         static_cast<ParsimonyScore::Weight>(0));
}

ParsimonyScore::Weight ParsimonyScore::AboveNode(
    ParsimonyScore::Weight edgeweight, ParsimonyScore::Weight childnodeweight) {
  return edgeweight + childnodeweight;
}
