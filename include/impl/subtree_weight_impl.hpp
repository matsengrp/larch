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
      dag_, dag_.GetDAG().GetRoot(), std::forward<WeightOps>(weight_ops),
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
  return SampleTreeImpl(
      std::forward<WeightOps>(weight_ops), [this, &weight_ops](auto clade) {
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
    optimum_edgeids.push_back(clade.at(i));
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

  ExtractTree(
      dag_, dag_.GetDAG().GetRoot(), std::forward<WeightOps>(weight_ops),
      [this, &distribution_maker](Node node, CladeIdx clade_idx) {
        auto clade = node.GetClade(clade_idx);
        Assert(not clade.empty());
        return clade.at(distribution_maker(clade)(random_generator_));
      },
      result, result_dag_ids);

  return {std::move(result), std::move(result_dag_ids)};
}

template <typename WeightOps>
template <typename EdgeSelector>
void SubtreeWeight<WeightOps>::ExtractTree(const MADAG& input_dag, Node node,
                                           WeightOps&& weight_ops,
                                           EdgeSelector&& edge_selector, MADAG& result,
                                           std::vector<NodeId>& result_dag_ids) {
  ComputeWeightBelow(node, std::forward<WeightOps>(weight_ops));
  CladeIdx clade_idx{0};

  NodeId parent_id{result.GetDAG().GetNodesCount()};
  result.AddNode(parent_id);
  GetOrInsert(result_dag_ids, parent_id) = node.GetId();

  if (not input_dag.GetCompactGenomes().empty()) {
    result.AppendCompactGenome(
        input_dag.GetCompactGenomes().at(node.GetId().value).Copy());
  }

  for (auto clade : node.GetClades()) {
    Edge edge = edge_selector(node, clade_idx);
    ++clade_idx.value;

    EdgeId edge_id{result.GetDAG().GetEdgesCount()};
    NodeId child_id{result.GetDAG().GetNodesCount()};

    result.AddEdge(edge_id, parent_id, child_id, edge.GetClade());

    if (not input_dag.GetEdgeMutations().empty()) {
      result.AppendEdgeMutations(
          input_dag.GetEdgeMutations().at(edge.GetId().value).Copy());
    }

    ExtractTree(dag_, edge.GetChild(), std::forward<WeightOps>(weight_ops),
                std::forward<EdgeSelector>(edge_selector), result, result_dag_ids);
  }

  if (node.IsRoot()) {
    result.BuildConnections();
  }

  for (auto node : result.GetDAG().GetNodes()) {
    size_t idx = node.GetId().value;
    const std::optional<std::string>& old_sample_id =
        input_dag.GetDAG().GetNodes().at(idx).GetSampleId();
    if (node.IsLeaf() and old_sample_id.has_value()) {
      node.SetSampleId(std::optional<std::string>{old_sample_id});
    }
  }
}
