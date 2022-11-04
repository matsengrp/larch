#include <algorithm>

template <typename BinaryOperatorWeightOps>
typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight
SimpleWeightOps<BinaryOperatorWeightOps>::ComputeLeaf(MADAG dag, NodeId node_id) {
  return binary_operator_weight_ops_.ComputeLeaf(dag, node_id);
}

template <typename BinaryOperatorWeightOps>
typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight
SimpleWeightOps<BinaryOperatorWeightOps>::ComputeEdge(MADAG dag, EdgeId edge_id) {
  return binary_operator_weight_ops_.ComputeEdge(dag, edge_id);
}

template <typename BinaryOperatorWeightOps>
std::pair<typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight,
          std::vector<size_t>>
SimpleWeightOps<BinaryOperatorWeightOps>::WithinCladeAccumOptimum(
    std::vector<typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight> inweights) {
  typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight optimal_weight =
      inweights[0];
  std::vector<size_t> optimal_indices;
  size_t inweight_idx = 0;
  for (auto weight : inweights) {
    if (binary_operator_weight_ops_.Compare(weight, optimal_weight)) {
      optimal_weight = weight;
      optimal_indices.clear();
      optimal_indices.push_back(inweight_idx);
    } else if (binary_operator_weight_ops_.CompareEqual(weight, optimal_weight)) {
      optimal_indices.push_back(inweight_idx);
    }
    inweight_idx++;
  }
  return {binary_operator_weight_ops_.Combine(BinaryOperatorWeightOps::Identity,
                                              optimal_weight),
          optimal_indices};
}

template <typename BinaryOperatorWeightOps>
typename BinaryOperatorWeightOps::Weight
SimpleWeightOps<BinaryOperatorWeightOps>::BetweenClades(
    std::vector<typename SimpleWeightOps<BinaryOperatorWeightOps>::Weight> inweights) {
  typename BinaryOperatorWeightOps::Weight result = BinaryOperatorWeightOps::Identity;
  for (auto weight : inweights) {
    result = binary_operator_weight_ops_.Combine(result, weight);
  }
  return result;
}

template <typename BinaryOperatorWeightOps>
typename BinaryOperatorWeightOps::Weight
SimpleWeightOps<BinaryOperatorWeightOps>::AboveNode(
    typename BinaryOperatorWeightOps::Weight edgeweight,
    typename BinaryOperatorWeightOps::Weight childnodeweight) {
  return binary_operator_weight_ops_.Combine(edgeweight, childnodeweight);
}
