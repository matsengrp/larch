/**
  Representation of edge mutations parsimony-based weight scoring for MADAG.

  This type is meant to be used as a parameter to SubtreeWeight.

 */

#pragma once

#include "mutation_annotated_dag.hpp"

struct NOP {
  template <typename... Args>
  void operator()(Args&&...) {}
};

template <typename MinWeightEdgeCallback = NOP, typename VisitNodeCallback = NOP>
struct ParsimonyScore {
  using Weight = size_t;
  constexpr static Weight MaxWeight = std::numeric_limits<size_t>::max();
  constexpr static Weight Identity = 0;
  Weight ComputeLeaf(Node node);
  Weight ComputeEdge(Edge edge);
  bool Compare(Weight lhs, Weight rhs);
  Weight Combine(Weight lhs, Weight rhs);

  // Callbacks
  void MinWeightEdge(Edge edge);
  void VisitNode(Node node, Weight&& weight);

  MADAG& dag_;
  MinWeightEdgeCallback min_weight_edge_callback_ = {};
  VisitNodeCallback visit_node_callback_ = {};
};

#include "impl/parsimony_score_impl.hpp"