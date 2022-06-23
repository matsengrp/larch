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
  Merge(std::string_view reference_sequence);

  Merge(Merge&&) = delete;
  Merge(const Merge&) = delete;
  Merge& operator=(Merge&&) = delete;
  Merge& operator=(const Merge&) = delete;

  void AddDAGs(const std::vector<std::reference_wrapper<MADAG>>& dags,
               bool have_compact_genomes = false);

  DAG& GetResult();
  const DAG& GetResult() const;
  const std::unordered_map<NodeLabel, NodeId>& GetResultNodes() const;
  [[nodiscard]] std::vector<EdgeMutations> ComputeResultEdgeMutations() const;

 private:
  void ComputeCompactGenomes(const std::vector<size_t>& tree_idxs);

  void ComputeLeafSets(const std::vector<size_t>& tree_idxs);

  void MergeTrees(const std::vector<size_t>& tree_idxs);

  static std::vector<LeafSet> ComputeLeafSets(const MADAG& dag,
                                              const std::vector<NodeLabel>& labels);

  std::string_view reference_sequence_;
  std::vector<std::reference_wrapper<MADAG>> trees_;

  ConcurrentUnorderedSet<CompactGenome> all_compact_genomes_;
  ConcurrentUnorderedSet<LeafSet> all_leaf_sets_;
  std::vector<std::vector<NodeLabel>> tree_labels_;

  std::unordered_map<NodeLabel, NodeId> result_nodes_;
  ConcurrentUnorderedMap<EdgeLabel, EdgeId> result_edges_;
  DAG result_dag_;
};
