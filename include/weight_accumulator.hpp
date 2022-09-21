/**
  Provides a class and templated WeightOps struct for counting the weights of all trees
  in a DAG.

  The WeightAccumulator struct is meant to be used as a parameter to SubtreeWeight.

  This should correctly accumulate weights for any WeightOps for which
  WeightOps::WithinCladeAccumOptimum called on a list containing a single Weight returns
  only that weight, and for which WeightOps::BetweenClades can be decomposed as a
  commutative binary operation on any list of Weights.

 */

#pragma once

#include "mutation_annotated_dag.hpp"
template <typename WeightOps>
class WeightCounter {
 public:
  WeightCounter(WeightCounter<WeightOps>&);
  WeightCounter(WeightOps&& weight_ops);
  WeightCounter(std::vector<typename WeightOps::Weight>, WeightOps&& weight_ops);
  WeightCounter(std::map<typename WeightOps::Weight, size_t>,
                const WeightOps&& weight_ops);
  /* A union of multisets */
  inline WeightCounter operator+(WeightCounter rhs);
  /* a cartesian product of multisets, applying
   * weight_ops.BetweenClades to pairs in the product */
  inline WeightCounter operator*(WeightCounter rhs);
  WeightCounter<WeightOps> operator=(WeightCounter<WeightOps> rhs);
  WeightCounter<WeightOps>& operator=(WeightCounter<WeightOps>& rhs);
  WeightCounter<WeightOps>&& operator=(WeightCounter<WeightOps>&& rhs);
  std::map<typename WeightOps::Weight, size_t> GetWeights();
  WeightOps&& GetWeightOps();

 private:
  std::map<typename WeightOps::Weight, size_t> weights_;
  const WeightOps&& weight_ops_;
};

template <typename WeightOps>
struct WeightAccumulator {
  using Weight = WeightCounter<WeightOps>;
  const WeightOps&& weight_ops_ = WeightOps();
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
  inline Weight AboveNode(Weight edgeweight, Weight childnodeweight);
};

#include "impl/weight_accumulator_impl.hpp"
