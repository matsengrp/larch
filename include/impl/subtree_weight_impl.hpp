#include <algorithm>

template <typename WeightOps>
SubtreeWeight<WeightOps>::SubtreeWeight(const DAG& dag)
    : weights_below_node_(dag.GetNodesCount(), WeightOps::Identity) {}

template <typename WeightOps>
typename WeightOps::Weight SubtreeWeight<WeightOps>::ComputeWeightBelow(
    Node node, WeightOps&& weight_ops) {
  typename WeightOps::Weight result = WeightOps::Identity;
  EdgeId min_weight_edge;
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
        min_weight_edge = edge;
      }
    }
    result = weight_ops.Combine(result, clade_min_weight);
    weight_ops.MinWeightEdge(node.GetDAG().Get(min_weight_edge));
  }
  cached = std::move(result);
done:
  weight_ops.VisitNode(node, std::forward<typename WeightOps::Weight>(cached));
  return cached;
}