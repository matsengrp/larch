#include <algorithm>

template <typename BinaryOperatorWeightOps>
typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight SimpleWeightOps<BinaryOperatorWeightOps>::ComputeLeaf(const MADAG& dag, NodeId node_id) {
  return BinaryOperatorWeightOps::ComputeLeaf(dag, node_id);
}

template <typename BinaryOperatorWeightOps>
typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight SimpleWeightOps<BinaryOperatorWeightOps>::ComputeEdge(const MADAG& dag, EdgeId edge_id) {
  return BinaryOperatorWeightOps::ComputeEdge(dag, edge_id);
}

template <typename BinaryOperatorWeightOps>
std::pair<typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight, std::vector<size_t>> SimpleWeightOps<BinaryOperatorWeightOps>::WithinCladeAccumOptimum(std::vector<typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight> inweights) {
    typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight optimal_weight = BinaryOperatorWeightOps::MaxWeight;
    std::vector<size_t> optimal_indices;
    size_t inweight_idx = 0;
    for (auto weight : inweights) {
        if (BinaryOperatorWeightOps::Compare(weight, optimal_weight)) {
            optimal_weight = weight;
            optimal_indices.clear();
            optimal_indices.push_back(inweight_idx);
        } else if (BinaryOperatorWeightOps::CompareEqual(weight, optimal_weight)) {
            optimal_indices.push_back(inweight_idx);
        }
        inweight_idx++;
    return {BinaryOperatorWeightOps::Combine(BinaryOperatorWeightOps::Identity, optimal_weight), optimal_indices};
    }
}

template <typename BinaryOperatorWeightOps>
typename BinaryOperatorWeightOps::Weight BetweenClades(std::vector<typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight> inweights) {
    typename BinaryOperatorWeightOps::Weight result = BinaryOperatorWeightOps::Identity;
    for (auto weight : inweights) {
        result = BinaryOperatorWeightOps::Combine(result, weight);
    }
    return result;
}

template <typename BinaryOperatorWeightOps>
typename BinaryOperatorWeightOps::Weight AboveNode(typename BinaryOperatorWeightOps::Weight edgeweight, typename BinaryOperatorWeightOps::Weight childnodeweight) {
    return BinaryOperatorWeightOps::Combine(edgeweight, childnodeweight);
}
