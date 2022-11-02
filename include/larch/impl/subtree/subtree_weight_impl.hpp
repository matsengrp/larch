#include <algorithm>
#include <type_traits>

template <typename WeightOps>
SubtreeWeight<WeightOps>::SubtreeWeight(const MADAG& dag)
    : dag_{dag},
      cached_weights_(dag_.GetDAG().GetNodesCount()),
      cached_min_weight_edges_(dag_.GetDAG().GetNodesCount()),
      random_device_{},
      random_generator_{random_device_()} {}

/*
 * Compute the Weight for the optimal subtree below `node`.
 *
 */
template <typename WeightOps>
typename WeightOps::Weight SubtreeWeight<WeightOps>::ComputeWeightBelow(
    Node node, WeightOps&& weight_ops) {
  NodeId node_id = node.GetId();
  auto& cached = cached_weights_.at(node_id.value);
  if (cached) {
    return *cached;
  }
  if (node.IsLeaf()) {
    cached = weight_ops.ComputeLeaf(dag_, node_id);
    return *cached;
  }
  std::vector<typename WeightOps::Weight> cladeweights;
  for (auto clade : node.GetClades()) {
    cladeweights.push_back(SubtreeWeight<WeightOps>::CladeWeight(
        clade, std::forward<WeightOps>(weight_ops)));
  }

  cached = weight_ops.BetweenClades(cladeweights);
  return *cached;
}

template <typename WeightOps>
MADAG SubtreeWeight<WeightOps>::TrimToMinWeight(WeightOps&& weight_ops) {
  MADAG result{dag_.GetReferenceSequence()};
  std::vector<NodeId> result_dag_ids;

  ExtractTree(
      dag_, dag_.GetDAG().GetRoot(), result.AppendNode(),
      std::forward<WeightOps>(weight_ops),
      [this](Node node, CladeIdx clade_idx) {
        return node.GetDAG().Get(
            // Probably don't want just the first one...?
            cached_min_weight_edges_.at(node.GetId().value).at(clade_idx.value)[0]);
      },
      result, result_dag_ids);

  return result;
}

template <typename WeightOps>
std::pair<MADAG, std::vector<NodeId>> SubtreeWeight<WeightOps>::SampleTree(
    WeightOps&& weight_ops) {
  return SampleTreeImpl(std::forward<WeightOps>(weight_ops), [](auto clade) {
    return std::uniform_int_distribution<size_t>{0, clade.size() - 1};
  });
}

struct TreeCount;
template <typename WeightOps>
std::pair<MADAG, std::vector<NodeId>> SubtreeWeight<WeightOps>::UniformSampleTree(
    WeightOps&& weight_ops) {
  static_assert(std::is_same_v<std::decay_t<WeightOps>, TreeCount>,
                "UniformSampleTree needs TreeCount");
  // Ensure cache is filled
  ComputeWeightBelow(dag_.GetDAG().GetRoot(), std::forward<WeightOps>(weight_ops));
  return SampleTreeImpl(std::forward<WeightOps>(weight_ops), [this](auto clade) {
    std::vector<double> probabilities;
    typename WeightOps::Weight sum{};
    for (NodeId child : clade | Transform::GetChild()) {
      sum += cached_weights_.at(child.value).value();
    }
    if (sum > 0) {
      for (NodeId child : clade | Transform::GetChild()) {
        probabilities.push_back(
            static_cast<double>(cached_weights_.at(child.value).value() / sum));
      }
    }
    return std::discrete_distribution<size_t>{probabilities.begin(),
                                              probabilities.end()};
  });
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
                             ComputeWeightBelow(dag_.GetDAG().Get(edge_id).GetChild(),
                                                std::forward<WeightOps>(weight_ops))));
  }
  auto clade_result = weight_ops.WithinCladeAccumOptimum(edge_weights);

  std::vector<EdgeId> optimum_edgeids;
  for (auto i : clade_result.second) {
    optimum_edgeids.push_back(
        clade.at(static_cast<ranges::range_difference_t<decltype(clade)>>(i)));
  }
  Edge first_edge = dag_.GetDAG().Get(clade[0]);
  GetOrInsert(GetOrInsert(cached_min_weight_edges_, first_edge.GetParentId().value),
              first_edge.GetClade().value) = optimum_edgeids;
  return clade_result.first;
}

template <typename WeightOps>
template <typename DistributionMaker>
std::pair<MADAG, std::vector<NodeId>> SubtreeWeight<WeightOps>::SampleTreeImpl(
    WeightOps&& weight_ops, DistributionMaker&& distribution_maker) {
  dag_.AssertUA();
  MADAG result{dag_.GetReferenceSequence()};
  std::vector<NodeId> result_dag_ids;

  NodeId parent_id = result.AppendNode();

  ExtractTree(
      dag_, dag_.GetDAG().GetRoot(), parent_id, std::forward<WeightOps>(weight_ops),
      [this, &distribution_maker](Node node, CladeIdx clade_idx) {
        auto clade = node.GetClade(clade_idx);
        Assert(not clade.empty());
        return clade.at(static_cast<ranges::range_difference_t<decltype(clade)>>(
            distribution_maker(clade)(random_generator_)));
      },
      result, result_dag_ids);

  result.BuildConnections();

  for (MutableNode node : result.GetDAG().GetNodes()) {
    const std::optional<std::string>& old_sample_id =
        dag_.GetDAG().Get(result_dag_ids.at(node.GetId().value)).GetSampleId();
    if (node.IsLeaf() and old_sample_id.has_value()) {
      node.SetSampleId(std::optional<std::string>{old_sample_id});
    }
  }

  return {std::move(result), std::move(result_dag_ids)};
}

template <typename WeightOps>
template <typename EdgeSelector>
void SubtreeWeight<WeightOps>::ExtractTree(const MADAG& input_dag, Node input_node,
                                           NodeId parent_id, WeightOps&& weight_ops,
                                           EdgeSelector&& edge_selector, MADAG& result,
                                           std::vector<NodeId>& result_dag_ids) {
  ComputeWeightBelow(input_node, std::forward<WeightOps>(weight_ops));
  CladeIdx clade_idx{0};

  GetOrInsert(result_dag_ids, parent_id) = input_node.GetId();

  if (not input_dag.GetCompactGenomes().empty()) {
    result.AppendCompactGenome(
        input_dag.GetCompactGenomes().at(input_node.GetId().value).Copy());
  }

  for (auto clade : input_node.GetClades()) {
    Assert(not clade.empty());

    Edge input_edge = edge_selector(input_node, clade_idx);
    ++clade_idx.value;

    NodeId child_id = result.AppendNode();
    result.AppendEdge(parent_id, child_id, input_edge.GetClade());

    if (not input_dag.GetEdgeMutations().empty()) {
      result.AppendEdgeMutations(
          input_dag.GetEdgeMutations().at(input_edge.GetId().value).Copy());
    }

    ExtractTree(dag_, input_edge.GetChild(), child_id,
                std::forward<WeightOps>(weight_ops),
                std::forward<EdgeSelector>(edge_selector), result, result_dag_ids);
  }
}
