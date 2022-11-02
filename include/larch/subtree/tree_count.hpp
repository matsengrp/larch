/**

  This type may be used as a parameter to SubtreeWeight for counting (sub)trees beneath
  each node.

 */

#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include "larch/mutation_annotated_dag.hpp"

struct TreeCount {
  using Weight = boost::multiprecision::cpp_int;
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

#include "larch/impl/subtree/tree_count_impl.hpp"
