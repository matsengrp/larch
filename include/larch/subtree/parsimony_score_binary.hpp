/**
  Representation of edge mutations parsimony-based weight scoring for MADAG.

  This type is meant to be used as a parameter to SimpleWeightOps, which can then be
  used as a parameter for SubtreeWeight.

 */

#pragma once

#include "larch/subtree/simple_weight_ops.hpp"

struct ParsimonyScore_ {
  using Weight = size_t;
  constexpr static Weight MaxWeight = std::numeric_limits<size_t>::max();
  constexpr static Weight Identity = 0;
  template <typename DAG>
  static Weight ComputeLeaf(DAG dag, NodeId node_id);
  template <typename DAG>
  static Weight ComputeEdge(DAG dag, EdgeId edge_id);
  static inline bool Compare(Weight lhs, Weight rhs);
  static inline bool CompareEqual(Weight lhs, Weight rhs);
  static inline bool IsIdentity(Weight weight);
  static inline Weight Combine(Weight lhs, Weight rhs);
};

struct MaxParsimonyScore_ : ParsimonyScore_ {
  constexpr static Weight MaxWeight = std::numeric_limits<size_t>::min();
  static inline bool Compare(Weight lhs, Weight rhs);
};

using BinaryParsimonyScore = SimpleWeightOps<ParsimonyScore_>;
using MaxBinaryParsimonyScore = SimpleWeightOps<MaxParsimonyScore_>;

#include "larch/impl/subtree/parsimony_score_binary_impl.hpp"
