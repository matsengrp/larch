/**
  Representation of edge mutations parsimony-based weight scoring for MADAG.

  This type is meant to be used as a parameter to SubtreeWeight.

 */

#pragma once

#include "larch/madag/mutation_annotated_dag.hpp"

struct ParsimonyScore {
  using Weight = size_t;

  template <typename DAG>
  static Weight ComputeLeaf(DAG dag, NodeId node_id);

  template <typename DAG>
  static Weight ComputeEdge(DAG dag, EdgeId edge_id);
  /*
   * Given a vector of weights for edges below a clade, compute the minimum
   * weight of them all, and return that minimum weight, and a vector
   * containing the indices of all elements of the passed vector that achieve
   * that minimum
   */
  inline static std::pair<Weight, std::vector<size_t>> WithinCladeAccumOptimum(
      const std::vector<Weight>&);
  /*
   * Given a vector of weights, one for each child clade, aggregate them
   */
  inline static Weight BetweenClades(const std::vector<Weight>&);
  inline static Weight AboveNode(Weight edgeweight, Weight childnodeweight);
};

#include "larch/impl/subtree/parsimony_score_impl.hpp"
