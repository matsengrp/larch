#pragma once

#include <string_view>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <shared_mutex>
#include <thread>
#include <atomic>

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for_each.h>

#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>

#include "mutation_annotated_dag.hpp"
#include "leaf_set.hpp"
#include "node_label.hpp"
#include "edge_label.hpp"

template <typename T>
using ConcurrentUnorderedSet =
    tbb::concurrent_unordered_set<T, std::hash<T>, std::equal_to<T>>;
template <typename K, typename V>
using ConcurrentUnorderedMap =
    tbb::concurrent_unordered_map<K, V, std::hash<K>, std::equal_to<K>>;

class Merge {
 public:
  /**
   * Construct a new Merge object, with the common reference sequence for all input
   * DAGs that will be merged later via the AddDAGs() method. The reference sequence is
   * externally owned, and should outlive the Merge object.
   */
  explicit Merge(std::string_view reference_sequence);

  Merge(Merge&&) = delete;
  Merge(const Merge&) = delete;
  Merge& operator=(Merge&&) = delete;
  Merge& operator=(const Merge&) = delete;

  /**
   * Add DAGs to be merged. The input DAGs are externally owned, and should outlive the
   * Merge object. If the have_compact_genomes parameter is false, the per-node compact
   * genomes of the input trees will be computed in parallel during the call to AddDAGs.
   * Otherwise the compact genomes stored in the DAGs will be used, and will be moved
   * into the Merge object storage to avoid duplication.
   */
  void AddDAGs(const std::vector<std::reference_wrapper<MADAG>>& dags,
               bool have_compact_genomes = false);

  /**
   * Get the DAG resulting from merge
   * @{
   */
  MADAG& GetResult();
  const MADAG& GetResult() const;
  /** @} */

  /**
   * Access the labels of the resulting DAG's nodes.
   */
  const std::unordered_map<NodeLabel, NodeId>& GetResultNodes() const;

  /**
   * Compute the mutations on the resulting DAG's edges. Can be used to build a MADAG
   * from the result.
   */
  [[nodiscard]] std::vector<EdgeMutations> ComputeResultEdgeMutations() const;

 private:
  void ComputeCompactGenomes(const std::vector<size_t>& tree_idxs);

  void ComputeLeafSets(const std::vector<size_t>& tree_idxs);

  void MergeTrees(const std::vector<size_t>& tree_idxs);

  static std::vector<LeafSet> ComputeLeafSets(const MADAG& dag,
                                              const std::vector<NodeLabel>& labels);

  // Vector of externally owned input DAGs.
  std::vector<std::reference_wrapper<MADAG>> trees_;

  // Every unique node compact genome, found among all input DAGs.
  ConcurrentUnorderedSet<CompactGenome> all_compact_genomes_;

  // Every unique node leaf set, found among all input DAGs.
  ConcurrentUnorderedSet<LeafSet> all_leaf_sets_;

  // Node labels for all input DAGs. Outer vector is indexed by input tree idx, inner
  // vector is indexed by node id.
  std::vector<std::vector<NodeLabel>> tree_labels_;

  // Node ids of the resulting DAG's nodes.
  std::unordered_map<NodeLabel, NodeId> result_nodes_;

  // Edge ids of the resulting DAG's edges.
  ConcurrentUnorderedMap<EdgeLabel, EdgeId> result_edges_;

  // Resulting DAG from merging the input DAGs.
  MADAG result_dag_;
};
