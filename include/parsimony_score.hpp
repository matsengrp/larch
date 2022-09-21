/**
  Representation of edge mutations parsimony-based weight scoring for MADAG.

  This type is meant to be used as a parameter to SubtreeWeight.

 */

#pragma once

#include "mutation_annotated_dag.hpp"

struct ParsimonyScore {
  using Weight = size_t;
  inline Weight ComputeLeaf(const MADAG& dag, NodeId node_id);
  inline Weight ComputeEdge(const MADAG& dag, EdgeId edge_id);
  /*
   * Given a vector of weights for edges below a clade, compute the minimum
   * weight of them all, and return that minimum weight, and a vector
   * containing the indices of all elements of the passed vector that achieve
   * that minimum
   */
  inline std::pair<Weight, std::vector<size_t>> WithinCladeAccumOptimum(
      std::vector<Weight>);
  /*
   * Given a vector of weights, one for each child clade, aggregate them
   */
  inline Weight BetweenClades(std::vector<Weight>);
  inline Weight AboveNode(Weight edgweight, Weight childnodeweight);
};

#include "impl/parsimony_score_impl.hpp"
