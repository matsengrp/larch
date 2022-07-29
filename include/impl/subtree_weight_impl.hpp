#include <algorithm>

template <typename T, typename WeightOps>
SubtreeWeight<T, WeightOps>::SubtreeWeight(const DAG& dag)
    : weights_below_node_(dag.GetNodesCount(), WeightOps::Identity) {}

template <typename T, typename WeightOps>
T SubtreeWeight<T, WeightOps>::ComputeWeightBelow(Node node, WeightOps&& weight_ops) {
  typename WeightOps::Weight result = WeightOps::Identity;
  auto& cached = weights_below_node_.at(node.GetId().value);
  if (cached != WeightOps::Identity) {
    goto done;
  }
  if (node.IsLeaf()) {
    cached = weight_ops.ComputeLeaf(node);
    goto done;
  }
  for (auto clade : node.GetClades()) {
    Assert(not clade.empty());
    typename WeightOps::Weight clade_min_weight = WeightOps::MaxWeight;
    for (auto edge : clade) {
      auto weight = weight_ops.Combine(
          weight_ops.ComputeEdge(edge),
          ComputeWeightBelow(edge.GetChild(), std::forward<WeightOps>(weight_ops)));
      if (weight_ops.Compare(weight, clade_min_weight)) {
        clade_min_weight = weight;
      }
    }
    result = weight_ops.Combine(result, clade_min_weight);
  }
  cached = std::move(result);
done:
  weight_ops.VisitNode(node, std::forward<typename WeightOps::Weight>(cached));
  return cached;
}