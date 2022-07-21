#include <algorithm>

template <typename T>
SubtreeWeight<T>::SubtreeWeight(const MADAG& dag) : dag_{dag} {
  weights_below_node_.resize(dag_.GetDAG().GetNodesCount());
}

template <typename T>
template <typename WeightOps>
T SubtreeWeight<T>::ComputeWeightBelow(Node node, WeightOps&& weight_ops) {
  if (node.IsLeaf()) {
    return weight_ops.ComputeLeaf(dag_, node);
  }
  typename WeightOps::Weight result = WeightOps::Identity;
  for (auto clade : node.GetClades()) {
    Assert(not clade.empty());
    typename WeightOps::Weight clade_min_weight = WeightOps::MaxWeight;
    for (auto edge : clade) {
      auto weight = weight_ops.Combine(
          weight_ops.ComputeEdge(dag_, edge),
          ComputeWeightBelow(edge.GetChild(), std::forward<WeightOps>(weight_ops)));
      if (weight_ops.Compare(weight, clade_min_weight)) {
        clade_min_weight = weight;
      }
    }
    result = weight_ops.Combine(result, clade_min_weight);
  }
  return result;
}