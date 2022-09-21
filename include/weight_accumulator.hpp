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
  WeightCounter(const WeightCounter<WeightOps>&);
  WeightCounter(WeightOps&& weight_ops);
  WeightCounter(const WeightOps& weight_ops) : weight_ops_{weight_ops} {}
  WeightCounter(const std::vector<typename WeightOps::Weight>& weights,
                const WeightOps& weight_ops);
  WeightCounter(const std::map<typename WeightOps::Weight, size_t>& weights,
                const WeightOps& weight_ops);
  /* A union of multisets */
  WeightCounter operator+(const WeightCounter& rhs) const;
  /* a cartesian product of multisets, applying
   * weight_ops.BetweenClades to pairs in the product */
  WeightCounter operator*(const WeightCounter& rhs) const;
  WeightCounter<WeightOps>& operator=(const WeightCounter<WeightOps>& rhs);
  WeightCounter<WeightOps>& operator=(WeightCounter<WeightOps>&& rhs);
  const std::map<typename WeightOps::Weight, size_t>& GetWeights() const;
  const WeightOps& GetWeightOps() const;

 private:
  std::map<typename WeightOps::Weight, size_t> weights_;
  WeightOps weight_ops_;
};

template <typename WeightOps>
struct WeightAccumulator {
  using Weight = WeightCounter<WeightOps>;
  WeightAccumulator() : weight_ops_{} {}
  WeightAccumulator(const WeightOps& ops) : weight_ops_{ops} {}
  WeightOps weight_ops_;
  inline Weight ComputeLeaf(const MADAG& dag, NodeId node_id);
  inline Weight ComputeEdge(const MADAG& dag, EdgeId edge_id);

  /*
   * Given a vector of weights for edges below a clade, compute the minimum
   * weight of them all, and return that minimum weight, and a vector
   * containing the indices of all elements of the passed vector that achieve
   * that minimum
   */
  inline std::pair<Weight, std::vector<size_t>> WithinCladeAccumOptimum(
      const std::vector<Weight>&);
  /*
   * Given a vector of weights, one for each child clade, aggregate them
   */
  inline Weight BetweenClades(const std::vector<Weight>&) const;
  inline Weight AboveNode(Weight edgeweight, Weight childnodeweight);
};

#include "impl/weight_accumulator_impl.hpp"
