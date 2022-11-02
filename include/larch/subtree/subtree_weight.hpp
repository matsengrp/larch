/**
  SubtreeWeight is used to calculate the (possibly optimal) tree weights achievable
below a given node in a DAG.

  SubtreeWeight expects as a template parameter a descriptor of required operations on
the weight type, including at least the values of this example struct:

struct ExampleWeightOps {
  using Weight = size_t;  // Provide any weight type
  inline Weight ComputeLeaf(const MADAG& dag, NodeId node_id);  // The value assigned to
each leaf node inline Weight ComputeEdge(const MADAG& dag, EdgeId edge_id);  // The
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

#include "larch/mutation_annotated_dag.hpp"

template <typename WeightOps>
class SubtreeWeight {
 public:
  explicit SubtreeWeight(const MADAG& dag);

  typename WeightOps::Weight ComputeWeightBelow(Node node, WeightOps&& weight_ops);

  [[nodiscard]] MADAG TrimToMinWeight(WeightOps&& weight_ops);

  [[nodiscard]] std::pair<MADAG, std::vector<NodeId>> SampleTree(
      WeightOps&& weight_ops);

  [[nodiscard]] std::pair<MADAG, std::vector<NodeId>> UniformSampleTree(
      WeightOps&& weight_ops);

 private:
  template <typename CladeRange>
  typename WeightOps::Weight CladeWeight(CladeRange&& clade, WeightOps&& weight_ops);

  template <typename DistributionMaker>
  [[nodiscard]] std::pair<MADAG, std::vector<NodeId>> SampleTreeImpl(
      WeightOps&& weight_ops, DistributionMaker&& distribution_maker);

  template <typename EdgeSelector>
  void ExtractTree(const MADAG& input_dag, Node node, NodeId parent_id,
                   WeightOps&& weight_ops, EdgeSelector&& edge_selector, MADAG& result,
                   std::vector<NodeId>& result_dag_ids);

  const MADAG& dag_;

  // Indexed by NodeId.
  std::vector<std::optional<typename WeightOps::Weight>> cached_weights_;

  // outermost vector indexed by NodeId, next vector indexed by CladeIdx, innermost
  // vector records which EdgeIds achieve minimum in that clade.
  std::vector<std::vector<std::vector<EdgeId>>> cached_min_weight_edges_;

  std::random_device random_device_;
  std::mt19937 random_generator_;
};

#include "larch/impl/subtree/subtree_weight_impl.hpp"
