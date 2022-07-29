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
  Weight ComputeLeaf(Node node);
  Weight ComputeEdge(Edge edge);
  bool Compare(Weight lhs, Weight rhs);
  Weight Combine(Weight lhs, Weight rhs);
  void VisitNode(Node node, Weight&& weight);

  MADAG& dag_;
};

#include "impl/parsimony_score_impl.hpp"