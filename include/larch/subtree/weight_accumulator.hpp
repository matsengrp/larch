/**
  Provides a class and templated WeightOps struct for counting the weights of all trees
  in a DAG.

  The WeightAccumulator struct is meant to be used as a parameter to SubtreeWeight.

  This should correctly accumulate weights for any WeightOps for which
  WeightOps::WithinCladeAccumOptimum called on a list containing a single Weight returns
  only that weight, and for which WeightOps::BetweenClades can be decomposed as a
  commutative binary operation on any list of Weights. (TODO this last requirement could
be avoided by another implementation)

  WeightAccumulator expects a WeightOps type with the same data/methods as the following
example:

  struct TreeWeightOps {
  using Weight = size_t;  // Provide any weight type
  inline Weight ComputeLeaf(const MADAG& dag, NodeId node_id);  // The value assigned to
each leaf node inline Weight ComputeEdge(const MADAG& dag, EdgeId edge_id);  // The
value assigned to each edge
  // Describes how to aggregate weights for alternative subtrees below a clade.
  // The returned pair contains the aggregated weight (for example, the optimal
  // one), and a vector containing indices for optimal weights in the input
  // vector. If optimality is undefined, all indices should be returned.
  inline std::pair<Weight, std::vector<size_t>>
  inline Weight BetweenClades(std::vector<Weight>);
  // Describes how to aggregate the weight of a subtree below a node, and the
  // weight of an edge above that node. The edge weight is always the first argument.
  inline Weight AboveNode(Weight edgweight, Weight childnodeweight);
};

  Note that unlike for the template parameter to SubtreeWeight, no
WithinCladeAccumOptimum method is needed, since each child clade in a tree may have only
one descendant edge.

 */

#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include "larch/mutation_annotated_dag.hpp"

using Count = boost::multiprecision::cpp_int;

template <typename WeightOps>
class WeightCounter {
 public:
  WeightCounter(const WeightCounter<WeightOps>&);
  WeightCounter(WeightOps&& weight_ops);
  WeightCounter(const WeightOps& weight_ops) : weight_ops_{weight_ops} {}
  WeightCounter(const std::vector<typename WeightOps::Weight>& weights,
                const WeightOps& weight_ops);
  WeightCounter(const std::map<typename WeightOps::Weight, Count>& weights,
                const WeightOps& weight_ops);
  /* A union of multisets */
  WeightCounter operator+(const WeightCounter& rhs) const;
  /* a cartesian product of multisets, applying
   * weight_ops.BetweenClades to pairs in the product */
  WeightCounter operator*(const WeightCounter& rhs) const;
  WeightCounter<WeightOps>& operator=(const WeightCounter<WeightOps>& rhs);
  WeightCounter<WeightOps>& operator=(WeightCounter<WeightOps>&& rhs);
  bool operator==(const WeightCounter<WeightOps>& rhs);
  bool operator!=(const WeightCounter<WeightOps>& rhs);
  const std::map<typename WeightOps::Weight, Count>& GetWeights() const;
  const WeightOps& GetWeightOps() const;

 private:
  std::map<typename WeightOps::Weight, Count> weights_;
  WeightOps weight_ops_;
};

template <typename WeightOps>
std::ostream& operator<<(std::ostream& os,
                         const WeightCounter<WeightOps>& weight_counter);

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

#include "larch/impl/subtree/weight_accumulator_impl.hpp"
