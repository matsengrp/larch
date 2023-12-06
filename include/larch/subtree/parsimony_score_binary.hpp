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
  Weight ComputeLeaf(DAG dag, NodeId node_id) const;
  template <typename DAG>
  Weight ComputeEdge(DAG dag, EdgeId edge_id) const;
  inline bool Compare(Weight lhs, Weight rhs) const;
  inline bool CompareEqual(Weight lhs, Weight rhs) const;
  inline bool IsIdentity(Weight weight) const;
  inline Weight Combine(Weight lhs, Weight rhs) const;
};

struct MaxParsimonyScore_ : ParsimonyScore_ {
  constexpr static Weight MaxWeight = std::numeric_limits<size_t>::min();
  inline bool Compare(Weight lhs, Weight rhs) const;
};

using BinaryParsimonyScore = SimpleWeightOps<ParsimonyScore_>;
using MaxBinaryParsimonyScore = SimpleWeightOps<MaxParsimonyScore_>;

#include "larch/impl/subtree/parsimony_score_binary_impl.hpp"
