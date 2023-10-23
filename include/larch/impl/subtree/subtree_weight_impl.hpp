#include <algorithm>
#include <type_traits>
#include <set>

struct SumRFDistance;

template <typename WeightOps, typename DAG>
SubtreeWeight<WeightOps, DAG>::SubtreeWeight(DAG dag)
    : dag_{dag},
      cached_weights_(dag_.GetNodesCount()),
      cached_subtree_counts_(dag_.GetNodesCount()),
      cached_min_weight_edges_(dag_.GetNodesCount()),
      random_generator_{random_device_()} {}

template <typename WeightOps, typename DAG>
DAG SubtreeWeight<WeightOps, DAG>::GetDAG() const {
  return dag_;
}

/*
 * Compute the Weight for the optimal subtree below `node`.
 *
 */
template <typename WeightOps, typename DAG>
typename WeightOps::Weight SubtreeWeight<WeightOps, DAG>::ComputeWeightBelow(
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
    cladeweights.push_back(CladeWeight(clade, std::forward<WeightOps>(weight_ops)));
  }

  cached = weight_ops.BetweenClades(cladeweights);
  return *cached;
}

template <typename WeightOps, typename DAG>
ArbitraryInt SubtreeWeight<WeightOps, DAG>::MinWeightCount(Node node,
                                                           WeightOps&& weight_ops) {
  NodeId node_id = node.GetId();
  auto& cached = cached_subtree_counts_.at(node_id.value);
  if (cached) {
    return *cached;
  }
  if (node.IsLeaf()) {
    cached = 1;
    return *cached;
  }
  ArbitraryInt nodecount = 1;
  ArbitraryInt cladecount = 0;
  // This populates cached_min_weight_edges_:
  SubtreeWeight<WeightOps, DAG>::ComputeWeightBelow(
      node, std::forward<WeightOps>(weight_ops));
  for (auto& clade : cached_min_weight_edges_.at(node_id.value)) {
    cladecount = 0;
    for (auto child_edge_id : clade) {
      Node child = dag_.Get(child_edge_id).GetChild();
      cladecount += SubtreeWeight<WeightOps, DAG>::MinWeightCount(
          child, std::forward<WeightOps>(weight_ops));
    }
    nodecount *= cladecount;
  }
  cached = nodecount;
  return *cached;
}

template <typename WeightOps, typename DAG>
typename SubtreeWeight<WeightOps, DAG>::Storage
SubtreeWeight<WeightOps, DAG>::TrimToMinWeight(WeightOps&& weight_ops) {
  Storage result{{}};
  result.View().SetReferenceSequence(dag_.GetReferenceSequence());
  ExtractTree(
      dag_.GetRoot(), result.View().AppendNode(), std::forward<WeightOps>(weight_ops),
      [this](Node node, CladeIdx clade_idx) {
        return node.GetDAG().Get(
            // Probably don't want just the first one...?
            cached_min_weight_edges_.at(node.GetId().value).at(clade_idx.value)[0]);
      },
      result.View());

  return result;
}

template <typename WeightOps, typename DAG>
typename SubtreeWeight<WeightOps, DAG>::SampledDAGStorage
SubtreeWeight<WeightOps, DAG>::SampleTree(WeightOps&& weight_ops,
                                          std::optional<NodeId> below) {
  Node below_node = below.has_value() ? dag_.Get(*below) : dag_.GetRoot();
  return SampleTreeImpl(
      std::forward<WeightOps>(weight_ops),
      [](auto clade) {
        return std::uniform_int_distribution<size_t>{0, clade.size() - 1};
      },
      below_node);
}

struct TreeCount;
template <typename WeightOps, typename DAG>
typename SubtreeWeight<WeightOps, DAG>::SampledDAGStorage
SubtreeWeight<WeightOps, DAG>::UniformSampleTree(WeightOps&& weight_ops,
                                                 std::optional<NodeId> below) {
  static_assert(std::is_same_v<std::decay_t<WeightOps>, TreeCount>,
                "UniformSampleTree needs TreeCount");
  Node below_node = below.has_value() ? dag_.Get(*below) : dag_.GetRoot();
  // Ensure cache is filled
  ComputeWeightBelow(below_node, std::forward<WeightOps>(weight_ops));
  return SampleTreeImpl(
      std::forward<WeightOps>(weight_ops),
      [this](auto clade) {
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
      },
      below_node);
}

template <typename WeightOps, typename DAG>
typename SubtreeWeight<WeightOps, DAG>::SampledDAGStorage
SubtreeWeight<WeightOps, DAG>::MinWeightSampleTree(WeightOps&& weight_ops,
                                                   std::optional<NodeId> below) {
  Node below_node = below.has_value() ? dag_.Get(*below) : dag_.GetRoot();
  // Ensure cache is filled
  ComputeWeightBelow(below_node, std::forward<WeightOps>(weight_ops));
  return SampleTreeImpl(
      std::forward<WeightOps>(weight_ops),
      [this](auto clade) {
        Edge first_edge = dag_.Get(clade.at(0));
        auto cached_clade = cached_min_weight_edges_.at(first_edge.GetParentId().value)
                                .at(first_edge.GetClade().value);
        std::set<EdgeId> min_weight_edges(cached_clade.begin(), cached_clade.end());
        std::vector<double> probabilities;
        for (EdgeId child_edge : clade) {
          if (min_weight_edges.count(child_edge)) {
            probabilities.push_back(1);
          } else {
            probabilities.push_back(0);
          }
        }
        return std::discrete_distribution<size_t>{probabilities.begin(),
                                                  probabilities.end()};
      },
      below_node);
}

template <typename WeightOps, typename DAG>
typename SubtreeWeight<WeightOps, DAG>::SampledDAGStorage
SubtreeWeight<WeightOps, DAG>::MinWeightUniformSampleTree(WeightOps&& weight_ops,
                                                          std::optional<NodeId> below) {
  Node below_node = below.has_value() ? dag_.Get(*below) : dag_.GetRoot();
  // Ensure cache is filled
  // (This also calls ComputeWeightBelow)
  ComputeWeightBelow(dag_.GetRoot(), std::forward<WeightOps>(weight_ops));
  return SampleTreeImpl(
      std::forward<WeightOps>(weight_ops),
      [this, &weight_ops](auto clade) {
        Edge first_edge = dag_.Get(clade.at(0));
        auto& cached_clade = cached_min_weight_edges_.at(first_edge.GetParentId().value)
                                 .at(first_edge.GetClade().value);
        std::set<EdgeId> min_weight_edges(cached_clade.begin(), cached_clade.end());
        std::vector<ArbitraryInt> min_weight_counts;
        ArbitraryInt sum = 0;
        for (EdgeId child_edge : clade) {
          if (min_weight_edges.count(child_edge)) {
            ArbitraryInt child_count = MinWeightCount(
                dag_.Get(child_edge).GetChild(), std::forward<WeightOps>(weight_ops));
            sum += child_count;
            min_weight_counts.push_back(child_count);
          } else {
            min_weight_counts.push_back(0);
          }
        }
        std::vector<double> probabilities;
        for (auto count : min_weight_counts) {
          probabilities.push_back(static_cast<double>(count / sum));
        }
        return std::discrete_distribution<size_t>{probabilities.begin(),
                                                  probabilities.end()};
      },
      below_node);
}
template <typename WeightOps, typename DAG>
template <typename CladeRange>
typename WeightOps::Weight SubtreeWeight<WeightOps, DAG>::CladeWeight(
    CladeRange&& clade, WeightOps&& weight_ops) {
  Assert(not clade.empty());
  std::vector<typename WeightOps::Weight> edge_weights;
  for (auto edge_id : clade) {
    edge_weights.push_back(
        weight_ops.AboveNode(weight_ops.ComputeEdge(dag_, edge_id),
                             ComputeWeightBelow(dag_.Get(edge_id).GetChild(),
                                                std::forward<WeightOps>(weight_ops))));
  }
  auto clade_result = weight_ops.WithinCladeAccumOptimum(edge_weights);

  std::vector<EdgeId> optimum_edgeids;
  optimum_edgeids.reserve(clade_result.second.size());
  for (auto i : clade_result.second) {
    optimum_edgeids.push_back(
        clade.at(static_cast<ranges::range_difference_t<decltype(clade)>>(i)));
  }
  Edge first_edge = dag_.Get(clade.at(0));
  GetOrInsert(GetOrInsert(cached_min_weight_edges_, first_edge.GetParentId().value),
              first_edge.GetClade().value) = optimum_edgeids;
  return clade_result.first;
}

template <typename WeightOps, typename DAG>
template <typename DistributionMaker>
typename SubtreeWeight<WeightOps, DAG>::SampledDAGStorage
SubtreeWeight<WeightOps, DAG>::SampleTreeImpl(WeightOps&& weight_ops,
                                              DistributionMaker&& distribution_maker,
                                              Node below) {
  Assert(not below.IsLeaf());
  dag_.AssertUA();
  SampledDAGStorage result{Storage{{}}};
  result.View().SetReferenceSequence(dag_.GetReferenceSequence());
  ExtractTree(
      below, result.View().AppendNode(), std::forward<WeightOps>(weight_ops),
      [this, &distribution_maker](Node node, CladeIdx clade_idx) {
        auto clade = node.GetClade(clade_idx);
        Assert(not clade.empty());
        return clade.at(static_cast<ranges::range_difference_t<decltype(clade)>>(
            distribution_maker(clade)(random_generator_)));
      },
      result.View());

  result.View().BuildConnections();

  for (auto node : result.View().GetNodes()) {
    const std::optional<std::string>& old_sample_id =
        dag_.Get(node.GetOriginalId()).GetSampleId();
    if (node.IsLeaf() and old_sample_id.has_value()) {
      node = SampleId{std::optional<std::string>{old_sample_id}};
    }
  }

  if (not below.IsUA()) {
    Assert(not below.GetCompactGenome().empty());
    EdgeMutations muts = CompactGenome::ToEdgeMutations(dag_.GetReferenceSequence(), {},
                                                        below.GetCompactGenome());
    result.View().AddUA(muts);
  }

  return result;
}

template <typename WeightOps, typename DAG>
template <typename NodeType, typename EdgeSelector, typename MutableDAGType>
void SubtreeWeight<WeightOps, DAG>::ExtractTree(NodeType input_node,
                                                NodeId result_node_id,
                                                WeightOps&& weight_ops,
                                                EdgeSelector&& edge_selector,
                                                MutableDAGType result) {
  ComputeWeightBelow(input_node, std::forward<WeightOps>(weight_ops));
  auto result_node = result.Get(result_node_id);
  if constexpr (decltype(result_node)::template contains_feature<MappedNodes>) {
    result_node.SetOriginalId(input_node.GetId());
  }
  result_node = input_node.GetCompactGenome().Copy();
  result_node = SampleId{input_node.GetSampleId()};

  CladeIdx clade_idx{0};
  for (auto clade : input_node.GetClades()) {
    Assert(not clade.empty());

    Edge input_edge = edge_selector(input_node, clade_idx);
    ++clade_idx.value;

    NodeId result_child_id = result.AppendNode();
    auto result_edge =
        result.AppendEdge(result_node_id, result_child_id, input_edge.GetClade());

    result_edge.SetEdgeMutations(input_edge.GetEdgeMutations().Copy());

    ExtractTree(input_edge.GetChild(), result_child_id,
                std::forward<WeightOps>(weight_ops),
                std::forward<EdgeSelector>(edge_selector), result);
  }
}
