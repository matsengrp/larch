/**
  Representation of edge mutations parsimony-based weight scoring for MADAG.

  This type is meant to be used as a parameter to SubtreeWeight.

 */

#pragma once

#include "mutation_annotated_dag.hpp"

struct ParsimonyScore {
  using Weight = size_t;
  constexpr static Weight MaxWeight = std::numeric_limits<size_t>::max();
  constexpr static Weight Identity = 0;
  inline Weight ComputeLeaf(const MADAG& dag, NodeId node_id);
  inline Weight ComputeEdge(const MADAG& dag, EdgeId edge_id);
  inline bool Compare(Weight lhs, Weight rhs);
  inline bool IsIdentity(Weight weight);
  inline Weight Combine(Weight lhs, Weight rhs);
};

#include "impl/parsimony_score_impl.hpp"
