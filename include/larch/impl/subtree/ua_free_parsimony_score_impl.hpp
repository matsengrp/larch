#include <algorithm>

template <typename DAG>
UAFreeParsimonyScore::Weight UAFreeParsimonyScore::ComputeLeaf(DAG, NodeId) {
  return 0;
}

template <typename DAG>
UAFreeParsimonyScore::Weight UAFreeParsimonyScore::ComputeEdge(DAG dag,
                                                               EdgeId edge_id) {
  if (dag.Get(edge_id).GetParent().IsUA()) {
    return 0;
  }
  return dag.Get(edge_id).GetEdgeMutations().size();
}

std::pair<UAFreeParsimonyScore::Weight, std::vector<size_t>>
UAFreeParsimonyScore::WithinCladeAccumOptimum(
    const std::vector<UAFreeParsimonyScore::Weight>& inweights) {
  UAFreeParsimonyScore::Weight optimal_weight = std::numeric_limits<size_t>::max();
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

UAFreeParsimonyScore::Weight UAFreeParsimonyScore::BetweenClades(
    const std::vector<UAFreeParsimonyScore::Weight>& inweights) {
  return std::accumulate(inweights.begin(), inweights.end(),
                         static_cast<UAFreeParsimonyScore::Weight>(0));
}

UAFreeParsimonyScore::Weight UAFreeParsimonyScore::AboveNode(
    UAFreeParsimonyScore::Weight edgeweight,
    UAFreeParsimonyScore::Weight childnodeweight) {
  return edgeweight + childnodeweight;
}
