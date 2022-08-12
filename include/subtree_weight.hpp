/**
  SubtreeWeight is used to calculate the minimum weight below a given node in a DAG.

  The two template parameters are respectively the type of the weight, and a
  descriptor of some arithmetic operations on the weight type. An example of such
  descriptor struct can be seen in parsimony_score.hpp.

 */

#pragma once

#include <functional>

#include "mutation_annotated_dag.hpp"

template <typename WeightOps>
class SubtreeWeight {
 public:
  explicit SubtreeWeight(const MADAG& dag);

  typename WeightOps::Weight ComputeWeightBelow(Node node, WeightOps&& weight_ops);

 private:
  const MADAG& dag_;
  std::vector<typename WeightOps::Weight> weights_below_node_;
};

#include "impl/subtree_weight_impl.hpp"