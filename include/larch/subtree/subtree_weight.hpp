/**
  SubtreeWeight is used to calculate the (possibly optimal) tree weights achievable
below a given node in a DAG.

  SubtreeWeight expects as a template parameter a descriptor of required operations on
the weight type, including at least the values of this example struct:

struct ExampleWeightOps {
  using Weight = size_t;  // Provide any weight type
  inline Weight ComputeLeaf(MADAG dag, NodeId node_id);  // The value assigned to
each leaf node inline Weight ComputeEdge(MADAG dag, EdgeId edge_id);  // The
value assigned to each edge
  // Describes how to aggregate weights for alternative subtrees below a clade.
  // The returned pair contains the aggregated weight (for example, the optimal
  // one), and a vector containing indices for optimal weights in the input
  // vector. If optimality is undefined, all indices should be returned.
  inline std::pair<Weight, std::vector<size_t>>
WithinCladeAccumOptimum(std::vector<Weight>);
  // Describes how to aggregate weights of subtrees below different child clades
  inline Weight BetweenClades(std::vector<Weight>);
  // Describes how to aggregate the weight of a subtree below a node, and the
  // weight of an edge above that node. The edge weight is always the first argument.
  inline Weight AboveNode(Weight edgweight, Weight childnodeweight);
};

 */

#pragma once

#include <functional>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>

#include "larch/madag/mutation_annotated_dag.hpp"

using ArbitraryInt = boost::multiprecision::cpp_int;

template <typename DAG>
struct SampledDAGStorage;

template <typename DAG>
using SampledDAGStorageBase =
    ExtendStorageType<SampledDAGStorage<DAG>, DAG, Extend::Nodes<MappedNodes>>;

template <typename DAG>
struct SampledDAGStorage : SampledDAGStorageBase<DAG> {
  static void Consume(DAG&&);
  static void FromView(const DAG&);
  static SampledDAGStorage EmptyDefault() {
    static_assert(DAG::role == Role::Storage);
    return SampledDAGStorage{DAG::EmptyDefault()};
  }

 private:
  friend SampledDAGStorageBase<DAG>;
  SampledDAGStorage(DAG&& target)
      : SampledDAGStorageBase<DAG>{std::forward<DAG>(target)} {}
};

template <typename WeightOps, typename DAG>
class SubtreeWeight {
 public:
  using Storage = std::remove_const_t<typename DAG::StorageType>;
  using SampledDAGStorage = ::SampledDAGStorage<Storage>;
  using MutableDAG = typename DAG::MutableType;
  using Node = typename DAG::NodeView;
  using Edge = typename DAG::EdgeView;
  explicit SubtreeWeight(DAG dag);

  DAG GetDAG() const;

  typename WeightOps::Weight ComputeWeightBelow(Node node, WeightOps&& weight_ops);

  ArbitraryInt MinWeightCount(Node node, WeightOps&& weight_ops);

  [[nodiscard]] Storage TrimToMinWeight(WeightOps&& weight_ops);

  [[nodiscard]] SampledDAGStorage SampleTree(
      WeightOps&& weight_ops, std::optional<NodeId> below = std::nullopt);

  [[nodiscard]] SampledDAGStorage UniformSampleTree(
      WeightOps&& weight_ops, std::optional<NodeId> below = std::nullopt);

  [[nodiscard]] SampledDAGStorage MinWeightSampleTree(
      WeightOps&& weight_ops, std::optional<NodeId> below = std::nullopt);

  [[nodiscard]] SampledDAGStorage MinWeightUniformSampleTree(
      WeightOps&& weight_ops, std::optional<NodeId> below = std::nullopt);

 private:
  template <typename CladeRange>
  typename WeightOps::Weight CladeWeight(CladeRange&& clade, WeightOps&& weight_ops);

  template <typename DistributionMaker>
  [[nodiscard]] SampledDAGStorage SampleTreeImpl(WeightOps&& weight_ops,
                                                 DistributionMaker&& distribution_maker,
                                                 Node below);

  template <typename NodeType, typename EdgeSelector, typename MutableDAGType>
  void ExtractTree(NodeType input_node, NodeId result_node_id, WeightOps&& weight_ops,
                   EdgeSelector&& edge_selector, MutableDAGType result);

  DAG dag_;

  // Indexed by NodeId.
  std::vector<std::optional<typename WeightOps::Weight>> cached_weights_;
  std::vector<std::optional<ArbitraryInt>> cached_subtree_counts_;

  // outermost vector indexed by NodeId, next vector indexed by CladeIdx, innermost
  // vector records which EdgeIds achieve minimum in that clade.
  std::vector<std::vector<std::vector<EdgeId>>> cached_min_weight_edges_;

  std::random_device random_device_;
  std::mt19937 random_generator_;
};

#include "larch/impl/subtree/subtree_weight_impl.hpp"
