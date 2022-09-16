/**
  A way to build a WeightOp using binary operations...

  This type is meant to be used as a parameter to SubtreeWeight.

 */

#pragma once

#include "mutation_annotated_dag.hpp"

template <typename BinaryOperatorWeightOps>
struct SimpleWeightOps {
    using Weight = typename BinaryOperatorWeightOps::Weight;
    BinaryOperatorWeightOps binary_operator_weight_ops;
    inline Weight ComputeLeaf(const MADAG& dag, NodeId node_id);
    inline Weight ComputeEdge(const MADAG& dag, EdgeId edge_id);

    /*
     * Given a vector of weights for edges below a clade, compute the minimum
     * weight of them all, and return that minimum weight, and a vector
     * containing the indices of all elements of the passed vector that achieve
     * that minimum
     */
    inline std::pair<Weight, std::vector<size_t>> WithinCladeAccumOptimum(std::vector<Weight>);
    /*
     * Given a vector of weights, one for each child clade, aggregate them
     */
    inline Weight BetweenClades(std::vector<Weight>);
    inline Weight AboveNode(std::vector<Weight>);
};

#include "impl/simple_weight_ops_impl.hpp"
