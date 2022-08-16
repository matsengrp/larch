#include <algorithm>

template <typename WeightOps>
SubtreeWeight<WeightOps>::SubtreeWeight(const MADAG& dag)
    : dag_{dag},
      weights_below_node_(dag_.GetDAG().GetNodesCount(), WeightOps::Identity) {}

template <typename WeightOps>
template <typename MinEdgeCallback>
typename WeightOps::Weight SubtreeWeight<WeightOps>::ComputeWeightBelow(
    Node node, WeightOps&& weight_ops, MinEdgeCallback&& min_edge_callback) {
  auto& cached = weights_below_node_.at(node.GetId().value);
  if (not weight_ops.IsIdentity(cached)) {
    return cached;
  }
  if (node.IsLeaf()) {
    cached = weight_ops.ComputeLeaf(dag_, node);
    return cached;
  }
  typename WeightOps::Weight result = WeightOps::Identity;
  EdgeId min_weight_edge;
  for (auto clade : node.GetClades()) {
    Assert(not clade.empty());
    typename WeightOps::Weight clade_min_weight = WeightOps::MaxWeight;
    for (auto edge : clade) {
      auto weight = weight_ops.Combine(
          weight_ops.ComputeEdge(dag_, edge),
          ComputeWeightBelow(edge.GetChild(), std::forward<WeightOps>(weight_ops),
                             std::forward<MinEdgeCallback>(min_edge_callback)));
      if (weight_ops.Compare(weight, clade_min_weight)) {
        clade_min_weight = weight;
        min_weight_edge = edge;
      }
    }
    result = weight_ops.Combine(result, clade_min_weight);
    if constexpr (not std::is_same_v<MinEdgeCallback, NOP>) {
      min_edge_callback(node.GetDAG().Get(min_weight_edge));
    }
  }
  cached = result;
  return cached;
}

template <typename WeightOps>
MADAG SubtreeWeight<WeightOps>::TrimToMinWeight(WeightOps&& weight_ops) {
  MADAG result;
  result.GetReferenceSequence() = dag_.GetReferenceSequence();

  auto trim_fn = [&](Edge edge) {
    result.GetDAG().AddNode(edge.GetParent());
    result.GetDAG().AddNode(edge.GetChild());
    result.GetDAG().AddEdge(edge, edge.GetParent(), edge.GetChild(), edge.GetClade());
  };

  ComputeWeightBelow(dag_.GetDAG().GetRoot(), std::forward<WeightOps>(weight_ops),
                     trim_fn);

  return result;
}