#pragma once

#include <string_view>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <shared_mutex>
#include <thread>
#include <atomic>
#include <numeric>

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for_each.h>

#include "larch/madag/mutation_annotated_dag.hpp"
#include "larch/merge/leaf_set.hpp"
#include "larch/merge/node_label.hpp"
#include "larch/merge/edge_label.hpp"

template <typename T>
using ConcurrentUnorderedSet =
    tbb::concurrent_unordered_set<T, std::hash<T>, std::equal_to<T>>;
template <typename K, typename V>
using ConcurrentUnorderedMap =
    tbb::concurrent_unordered_map<K, V, std::hash<K>, std::equal_to<K>>;

class Merge {
 public:
  using Node = MADAG::Node;
  using Edge = MADAG::Edge;
  /**
   * Construct a new Merge object, with the common reference sequence for all input
   * DAGs that will be merged later via the AddDAGs() method. The reference sequence is
   * externally owned, and should outlive the Merge object.
   */
  inline explicit Merge(std::string_view reference_sequence);

  Merge(Merge&&) = delete;
  Merge(const Merge&) = delete;
  Merge& operator=(Merge&&) = delete;
  Merge& operator=(const Merge&) = delete;
  ~Merge() = default;

  /**
   * Add DAGs to be merged. The input DAGs are externally owned, and should outlive the
   * Merge object. If the have_compact_genomes parameter is false, the per-node compact
   * genomes of the input trees will be computed in parallel during the call to AddDAGs.
   * Otherwise the compact genomes stored in the DAGs will be used, and will be moved
   * into the Merge object storage to avoid duplication.
   */
  inline void AddDAGs(const std::vector<MADAG>& dags);

  /**
   * Get the DAG resulting from merge
   * @{
   */
  inline MADAG GetResult() const;
  /** @} */

  /**
   * Access the labels of the resulting DAG's nodes.
   */
  inline const std::unordered_map<NodeLabel, NodeId>& GetResultNodes() const;

  inline const std::vector<NodeLabel>& GetResultNodeLabels() const;

  /**
   * Compute the mutations on the resulting DAG's edges and store in the result MADAG.
   */
  inline void ComputeResultEdgeMutations();

  inline bool ContainsLeafset(const LeafSet& leafset) const;

 private:
  inline void ComputeLeafSets(const std::vector<size_t>& tree_idxs);

  inline void MergeTrees(const std::vector<size_t>& tree_idxs);

  inline static std::vector<LeafSet> ComputeLeafSets(
      MADAG dag, const std::vector<NodeLabel>& labels);

  // Vector of externally owned input DAGs.
  std::vector<MADAG> trees_;

  // Every unique node compact genome, found among all input DAGs.
  ConcurrentUnorderedSet<CompactGenome> all_compact_genomes_;

  // Every unique node leaf set, found among all input DAGs.
  ConcurrentUnorderedSet<LeafSet> all_leaf_sets_;

  // Node labels for all input DAGs. Outer vector is indexed by input tree idx, inner
  // vector is indexed by node id.
  std::vector<std::vector<NodeLabel>> tree_labels_;

  // Node ids of the resulting DAG's nodes.
  std::unordered_map<NodeLabel, NodeId> result_nodes_;
  std::vector<NodeLabel> result_node_labels_;

  // Edge ids of the resulting DAG's edges.
  ConcurrentUnorderedMap<EdgeLabel, EdgeId> result_edges_;

  // Resulting DAG from merging the input DAGs.
  MADAGStorage result_dag_storage_;
  MutableMADAG result_dag_;
};

#include "larch/impl/merge/merge_impl.hpp"