#include <algorithm>

template <typename WeightOps>
SubtreeWeight<WeightOps>::SubtreeWeight(const MADAG& dag)
    : dag_{dag},
      weights_below_node_(dag_.GetDAG().GetNodesCount(), WeightOps::Identity) {}

template <typename WeightOps>
typename WeightOps::Weight SubtreeWeight<WeightOps>::ComputeWeightBelow(
    Node node, WeightOps&& weight_ops) {
  auto& cached = weights_below_node_.at(node.GetId().value);
  if (cached != WeightOps::Identity) {
    return cached;
  }
  if (node.IsLeaf()) {
    cached = weight_ops.ComputeLeaf(dag_, node);
    return cached;
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
  cached = result;
  return cached;
}