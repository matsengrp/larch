#include <algorithm>

template <typename WeightOps>
SubtreeWeight<WeightOps>::SubtreeWeight(const MADAG& dag)
    : dag_{dag},
      cached_weights_(dag_.GetDAG().GetNodesCount()),
      weight_is_cached_(dag_.GetDAG().GetNodesCount(), false),
      random_device_{},
      random_generator_{random_device_()} {}

/*
 * Compute the Weight for the optimal subtree below `node`.
 *
 */
template <typename WeightOps>
typename WeightOps::Weight SubtreeWeight<WeightOps>::ComputeWeightBelow(Node node, WeightOps&& weight_ops) {
  NodeId node_id = node.GetId();
  auto& cached = cached_weights_.at(node_id.value);
  bool is_cached = weight_is_cached_.at(node_id.value);
  if (is_cached) {
    return cached;
  }
  if (node.IsLeaf()) {
    cached = weight_ops.ComputeLeaf(dag_, node_id);
    weight_is_cached_.at(node_id.value) = true;
    return cached;
  }
  std::vector<typename WeightOps::Weight> cladeweights;
  for (auto clade : node.GetClades()) {
      cladeweights.push_back(SubtreeWeight<WeightOps>::CladeWeight(clade, std::forward<WeightOps>(weight_ops)));
  }
  typename WeightOps::Weight result = weight_ops.BetweenClades(cladeweights);

  cached = result;
  weight_is_cached_.at(node_id.value) = true;
  return cached;
}

template <typename WeightOps>
MADAG SubtreeWeight<WeightOps>::TrimToMinWeight(WeightOps&& weight_ops) {
  MADAG result;
  result.GetReferenceSequence() = dag_.GetReferenceSequence();

  ExtractTree(
      dag_, dag_.GetDAG().GetRoot(), std::forward<WeightOps>(weight_ops),
      [this](Node node, CladeIdx clade_idx) {
        return node.GetDAG().Get(
            // Probably don't want just the first one...?
            cached_min_weight_edges_.at(node.GetId().value).at(clade_idx.value)[0]);
      },
      result);

  return result;
}

template <typename WeightOps>
MADAG SubtreeWeight<WeightOps>::SampleTree(WeightOps&& weight_ops) {
  MADAG result;
  result.GetReferenceSequence() = dag_.GetReferenceSequence();

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
typename WeightOps::Weight SubtreeWeight<WeightOps>::CladeWeight(
    CladeRange&& clade, WeightOps&& weight_ops) {
  assert(not clade.empty());
  std::vector<typename WeightOps::Weight> edge_weights;
  for (auto edge_id : clade) {
      edge_weights.push_back(
          weight_ops.AboveNode(weight_ops.ComputeEdge(dag_, edge_id),
          ComputeWeightBelow(dag_.GetDAG().Get(edge_id).GetChild(), std::forward<WeightOps>(weight_ops)))
      );
  }
  auto clade_result = weight_ops.WithinCladeAccumOptimum(edge_weights);

  std::vector<EdgeId> optimum_edgeids;
  for (auto i : clade_result.second) {
      optimum_edgeids.push_back(clade.at(i));
  }
  Edge first_edge = dag_.GetDAG().Get(clade[0]);
  cached_min_weight_edges_.at(first_edge.GetParentId().value).at(first_edge.GetClade().value) = optimum_edgeids;
  return clade_result.first;
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
  result.GetDAG().AddNode(parent_id);

  if (not input_dag.GetCompactGenomes().empty()) {
    result.GetCompactGenomes().push_back(
        input_dag.GetCompactGenomes().at(node.GetId().value).Copy());
  }

  for (auto clade : node.GetClades()) {
    Edge edge = edge_selector(node, clade_idx);
    ++clade_idx.value;

    EdgeId edge_id{result.GetDAG().GetEdgesCount()};
    NodeId child_id{result.GetDAG().GetNodesCount()};

    result.GetDAG().AddEdge(edge_id, parent_id, child_id, edge.GetClade());

    if (not input_dag.GetEdgeMutations().empty()) {
      result.GetEdgeMutations().push_back(
          input_dag.GetEdgeMutations().at(edge.GetId().value).Copy());
    }

    ExtractTree(dag_, edge.GetChild(), std::forward<WeightOps>(weight_ops),
                std::forward<EdgeSelector>(edge_selector), result);
  }

  if (node.IsRoot()) {
    result.GetDAG().BuildConnections();
  }
  for (auto node : input_dag.GetDAG().GetNodes()) {
      if (node.IsLeaf()) {
          if(node.GetSampleId()){
          std::cout << *node.GetSampleId() << "\n";
          }
      }
  }
  for (auto node : result.GetDAG().GetNodes()) {
      size_t idx = node.GetId().value;
      std::optional<std::string> old_sample_id = input_dag.GetDAG().GetNodes().at(idx).GetSampleId();
      if (node.IsLeaf() and (bool) old_sample_id) {
          std::cout << "nodeId: " << std::to_string(idx) << old_sample_id.value() << "\n";
          node.SetSampleId(old_sample_id);
      }
  }
}
