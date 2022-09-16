/**
  Representation of edge mutations parsimony-based weight scoring for MADAG.

  This type is meant to be used as a parameter to SubtreeWeight.

 */

#pragma once

#include "mutation_annotated_dag.hpp"
template <typename WeightOps>
class WeightCounter {
    public:
        WeightCounter();
        WeightCounter(std::vector<typename WeightOps::Weight>);
        WeightCounter operator+(WeightCounter lhs, WeightCounter rhs);
        WeightCounter operator*(WeightCounter lhs, WeightCounter rhs);
    private:
        std::map<typename WeightOps::Weight, size_t> weights_;
};

template <typename WeightOps>
struct WeightAccumulator {
  using Weight = WeightCounter<WeightOps>;
  /* constexpr static Weight MaxWeight = std::numeric_limits<size_t>::max(); */ 
  /* constexpr static Weight Identity = 0; */
  inline Weight ComputeLeaf(const MADAG& dag, NodeId node_id);
  inline Weight ComputeEdge(const MADAG& dag, EdgeId edge_id);
  /* inline bool Compare(Weight lhs, Weight rhs); */
  /* inline bool IsIdentity(Weight weight); */
  /* inline Weight Combine(Weight lhs, Weight rhs); */

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
  inline Weight AboveNode(Weight edgeweight, Weight childnodeweight);
};


#include "impl/parsimony_score_impl.hpp"
