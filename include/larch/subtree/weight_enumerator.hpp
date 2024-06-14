#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include "larch/madag/mutation_annotated_dag.hpp"

using Count = boost::multiprecision::cpp_int;

// template <typename WeightOps>
// struct WeightEnumerator {
//   using Weight = WeightCounter<WeightOps>;
//   WeightAccumulator() = default;
//   WeightAccumulator(const WeightOps& ops);

//   template <typename DAG>
//   Weight ComputeLeaf(DAG dag, NodeId node_id) const;

//   template <typename DAG>
//   Weight ComputeEdge(DAG dag, EdgeId edge_id) const;

//   /*
//    * Given a vector of weights for edges below a clade, compute the minimum
//    * weight of them all, and return that minimum weight, and a vector
//    * containing the indices of all elements of the passed vector that achieve
//    * that minimum
//    */
//   inline std::pair<Weight, std::vector<size_t>> WithinCladeAccumOptimum(
//       const std::vector<Weight>&) const;
//   /*
//    * Given a vector of weights, one for each child clade, aggregate them
//    */
//   inline Weight BetweenClades(const std::vector<Weight>&) const;
//   inline Weight AboveNode(Weight edgeweight, Weight childnodeweight) const;

//  private:
//   WeightOps weight_ops_ = {};
// };

#include "larch/impl/subtree/weight_enumerator_impl.hpp"
