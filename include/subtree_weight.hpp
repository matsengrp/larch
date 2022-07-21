#pragma once

#include <functional>

#include "mutation_annotated_dag.hpp"

template <typename T>
class SubtreeWeight {
 public:
  explicit SubtreeWeight(const MADAG& dag);

  template <typename WeightOps>
  T ComputeWeightBelow(Node node, WeightOps&& weight_ops);

 private:
  const MADAG& dag_;
  std::vector<T> weights_below_node_;
};

#include "impl/subtree_weight_impl.hpp"