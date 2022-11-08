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
  inline Weight ComputeLeaf(MADAG dag, NodeId node_id);
  inline Weight ComputeEdge(MADAG dag, EdgeId edge_id);
  inline bool Compare(Weight lhs, Weight rhs);
  inline bool CompareEqual(Weight lhs, Weight rhs);
  inline bool IsIdentity(Weight weight);
  inline Weight Combine(Weight lhs, Weight rhs);
};

using ParsimonyScore = SimpleWeightOps<ParsimonyScore_>;

#include "larch/impl/subtree/parsimony_score_binary_impl.hpp"
