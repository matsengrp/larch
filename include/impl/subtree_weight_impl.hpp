#include <algorithm>

template <typename WeightOps>
SubtreeWeight<WeightOps>::SubtreeWeight(const MADAG& dag)
    : dag_{dag},
      cached_weights_(dag_.GetDAG().GetNodesCount(), WeightOps::Identity),
      random_device_{},
      random_generator_{random_device_()} {}

template <typename WeightOps>
typename WeightOps::Weight SubtreeWeight<WeightOps>::ComputeWeightBelow(
    Node node, WeightOps&& weight_ops) {
  auto& cached = cached_weights_.at(node.GetId().value);
  if (not weight_ops.IsIdentity(cached)) {
    return cached;
  }
  if (node.IsLeaf()) {
    cached = weight_ops.ComputeLeaf(dag_, node);
    return cached;
  }

  typename WeightOps::Weight result = WeightOps::Identity;
  CladeIdx clade_idx{0};
  for (auto clade : node.GetClades()) {
    Assert(not clade.empty());
    auto min_weight = MinCladeWeight(clade, std::forward<WeightOps>(weight_ops));
    result = weight_ops.Combine(result, min_weight.first);
    auto& min_edges = GetOrInsert(cached_min_weight_edges_, clade_idx.value++);
    if (min_edges.empty()) {
      min_edges.resize(dag_.GetDAG().GetNodesCount());
    }
    min_edges.at(node.GetId().value) = min_weight.second;
  }

  cached = result;
  return cached;
}

template <typename WeightOps>
MADAG SubtreeWeight<WeightOps>::TrimToMinWeight(WeightOps&& weight_ops) {
  MADAG result{dag_.GetReferenceSequence()};

  ExtractTree(
      dag_, dag_.GetDAG().GetRoot(), std::forward<WeightOps>(weight_ops),
      [this](Node node, CladeIdx clade_idx) {
        return node.GetDAG().Get(
            cached_min_weight_edges_.at(clade_idx.value).at(node.GetId().value));
      },
      result);

  return result;
}

template <typename WeightOps>
MADAG SubtreeWeight<WeightOps>::SampleTree(WeightOps&& weight_ops) {
  MADAG result{dag_.GetReferenceSequence()};

  ExtractTree(
      dag_, dag_.GetDAG().GetRoot(), std::forward<WeightOps>(weight_ops),
      [this](Node node, CladeIdx clade_idx) {
        auto clade = node.GetClade(clade_idx);
        Assert(not clade.empty());
        std::uniform_int_distribution<size_t> distribuition{0, clade.size() - 1};
        return clade.at(distribuition(random_generator_));
      },
      result);

  return result;
}

template <typename WeightOps>
template <typename CladeRange>
std::pair<typename WeightOps::Weight, EdgeId> SubtreeWeight<WeightOps>::MinCladeWeight(
    CladeRange&& clade, WeightOps&& weight_ops) {
  typename WeightOps::Weight result = WeightOps::Identity;
  EdgeId min_weight_edge;
  Assert(not clade.empty());
  typename WeightOps::Weight clade_min_weight = WeightOps::MaxWeight;
  for (auto edge : clade) {
    auto weight = weight_ops.Combine(
        weight_ops.ComputeEdge(dag_, edge),
        ComputeWeightBelow(edge.GetChild(), std::forward<WeightOps>(weight_ops)));
    if (weight_ops.Compare(weight, clade_min_weight)) {
      clade_min_weight = weight;
      min_weight_edge = edge;
    }
  }
  result = weight_ops.Combine(result, clade_min_weight);

  return {result, min_weight_edge};
}

template <typename WeightOps>
template <typename EdgeSelector>
void SubtreeWeight<WeightOps>::ExtractTree(const MADAG& input_dag, Node node,
                                           WeightOps&& weight_ops,
                                           EdgeSelector&& edge_selector,
                                           MADAG& result) {
  ComputeWeightBelow(node, std::forward<WeightOps>(weight_ops));
  CladeIdx clade_idx{0};

  NodeId parent_id{result.GetDAG().GetNodesCount()};
  result.AddNode(parent_id);

  if (not input_dag.GetCompactGenomes().empty()) {
    result.AppendCompactGenome(
        input_dag.GetCompactGenomes().at(node.GetId().value).Copy());
  }

  for (auto clade : node.GetClades()) {
    Edge edge = edge_selector(node, clade_idx);
    clade_idx.value++;

    EdgeId edge_id{result.GetDAG().GetEdgesCount()};
    NodeId child_id{result.GetDAG().GetNodesCount()};

    result.AddEdge(edge_id, parent_id, child_id, edge.GetClade());

    if (not input_dag.GetEdgeMutations().empty()) {
      result.AppendEdgeMutations(
          input_dag.GetEdgeMutations().at(edge.GetId().value).Copy());
    }

    ExtractTree(dag_, edge.GetChild(), std::forward<WeightOps>(weight_ops),
                std::forward<EdgeSelector>(edge_selector), result);
  }

  if (node.IsRoot()) {
    result.BuildConnections();
  }
}
