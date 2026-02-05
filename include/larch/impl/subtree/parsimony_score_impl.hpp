#include <algorithm>

template <typename DAG>
ParsimonyScore::Weight ParsimonyScore::ComputeLeaf(DAG, NodeId) {
  return 0;
}

template <typename DAG>
ParsimonyScore::Weight ParsimonyScore::ComputeEdge(DAG dag, EdgeId edge_id) {
  if (dag.Get(edge_id).GetChild().IsLeaf()) {
    size_t mut_count = 0;
    for (auto mut : dag.Get(edge_id).GetEdgeMutations()) {
      if (not mut.second.second.IsAmbiguous() or mut.second.first.IsAmbiguous()) {
        mut_count += 1;
      }
    }
    return mut_count;
  }
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
