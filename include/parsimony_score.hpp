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
  Weight ComputeLeaf(const MADAG& dag, NodeId node_id);
  Weight ComputeEdge(const MADAG& dag, EdgeId edge_id);
  bool Compare(Weight lhs, Weight rhs);
  Weight Combine(Weight lhs, Weight rhs);
};

#include <algorithm>

ParsimonyScore::Weight ParsimonyScore::ComputeLeaf(const MADAG& dag, NodeId node_id) {
  return 0;
}

ParsimonyScore::Weight ParsimonyScore::ComputeEdge(const MADAG& dag, EdgeId edge_id) {
  return dag.GetEdgeMutations(edge_id).size();
}

bool ParsimonyScore::Compare(Weight lhs, Weight rhs) {
  return lhs < rhs;
}

ParsimonyScore::Weight ParsimonyScore::Combine(Weight lhs, Weight rhs) {
  return lhs + rhs;
}