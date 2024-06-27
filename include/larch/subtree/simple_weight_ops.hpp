/**
  A way to build a WeightOp type possessing all the values required by SubtreeWeight,
using simpler-to-define binary operations.

  The template parameter BinaryOperatorWeightOps is expected to provide at least the
values in the following example struct:

struct ExampleBinaryOperatorWeightOps {
  // Provide any weight type you like
  using Weight = size_t;
  // The identity element under Combine
  constexpr static Weight Identity = 0;
  // The value assigned to a leaf
  inline Weight ComputeLeaf(MADAG dag, NodeId node_id);
  // The value assigned to an edge
  node inline Weight ComputeEdge(MADAG dag, EdgeId edge_id);
  // The ordering operator: is lhs 'better than' rhs
  inline bool Compare(Weight lhs, Weight rhs);
  // A custom implemented equality
  inline bool CompareEqual(Weight lhs, Weight rhs);
  // A binary operation that respects ordering, e.g. '+'
  inline Weight Combine(Weight lhs, Weight rhs);
};

 */

#pragma once

#include "larch/madag/mutation_annotated_dag.hpp"

template <typename BinaryOperatorWeightOps>
struct SimpleWeightOps {
  using Weight = typename BinaryOperatorWeightOps::Weight;

  SimpleWeightOps() = default;

  explicit SimpleWeightOps(BinaryOperatorWeightOps&& ops)
      : binary_operator_weight_ops_{std::forward<BinaryOperatorWeightOps>(ops)} {}

  template <typename DAG>
  Weight ComputeLeaf(DAG dag, NodeId node_id) const;

  template <typename DAG>
  Weight ComputeEdge(DAG dag, EdgeId edge_id) const;
  inline std::pair<Weight, std::vector<size_t>> WithinCladeAccumOptimum(
      std::vector<Weight>) const;
  inline Weight BetweenClades(std::vector<Weight>) const;
  inline Weight AboveNode(Weight edgeweight, Weight childnodeweight) const;

  inline const BinaryOperatorWeightOps& GetOps() const;

 private:
  BinaryOperatorWeightOps binary_operator_weight_ops_ = BinaryOperatorWeightOps();
};

#include "larch/impl/subtree/simple_weight_ops_impl.hpp"
