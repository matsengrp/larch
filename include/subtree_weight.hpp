/**
  SubtreeWeight is used to calculate the minimum weight below a given node in a DAG.

  The two template parameters are respectively the type of the weight, and a
  descriptor of some arithmetic operations on the weight type. An example of such
  descriptor struct can be seen in parsimony_score.hpp.

 */

#pragma once

#include <functional>
#include <random>

#include "mutation_annotated_dag.hpp"

template <typename WeightOps>
class SubtreeWeight {
 public:
  explicit SubtreeWeight(const MADAG& dag);

  typename WeightOps::Weight ComputeWeightBelow(Node node, WeightOps&& weight_ops);

  [[nodiscard]] MADAG TrimToMinWeight(WeightOps&& weight_ops);

  [[nodiscard]] MADAG SampleTree(WeightOps&& weight_ops);

 private:
  template <typename CladeRange>
  std::pair<typename WeightOps::Weight, EdgeId> MinCladeWeight(CladeRange&& clade,
                                                               WeightOps&& weight_ops);

  template <typename EdgeSelector>
  void ExtractTree(const MADAG& input_dag, Node node, WeightOps&& weight_ops,
                   EdgeSelector&& edge_selector, MADAG& result);

  const MADAG& dag_;

  // Indexed by NodeId.
  std::vector<typename WeightOps::Weight> cached_weights_;

  // Outer vector indexed by CladeIdx, inner vector indexed by NodeId.
  std::vector<std::vector<EdgeId>> cached_min_weight_edges_;

  std::random_device random_device_;
  std::mt19937 random_generator_;
};

#include "impl/subtree_weight_impl.hpp"