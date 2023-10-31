#pragma once

#include <string_view>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <shared_mutex>
#include <mutex>
#include <thread>
#include <atomic>
#include <numeric>

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_vector.h>

#include "larch/madag/mutation_annotated_dag.hpp"
#include "larch/merge/leaf_set.hpp"
#include "larch/merge/node_label.hpp"
#include "larch/merge/edge_label.hpp"
#include "larch/parallel/parallel_common.hpp"

template <typename T>
using ConcurrentUnorderedSet =
    tbb::concurrent_unordered_set<T, std::hash<T>, std::equal_to<T>>;

template <typename K, typename V>
using ConcurrentUnorderedMap = SharedState<std::unordered_map<K, V>>;

template <typename Target = DefaultDAGStorage>
struct MergeDAGStorage;

template <typename Target>
struct LongNameOf<MergeDAGStorage<Target>> {
  using type = ExtendStorageType<
      MergeDAGStorage<Target>, DefaultDAGStorage,
      Extend::Nodes<Deduplicate<CompactGenome>, Deduplicate<SampleId>>,
      Extend::Edges<EdgeMutations>, Extend::DAG<ReferenceSequence>>;
};

template <typename Target>
struct MergeDAGStorage : LongNameOf<MergeDAGStorage<Target>>::type {
  SHORT_NAME(MergeDAGStorage);
};

using MergeDAG = DAGView<const MergeDAGStorage<>>;
using MutableMergeDAG = DAGView<MergeDAGStorage<>>;

class Merge {
  using NodeLabelsContainer = ContiguousMap<NodeId, NodeLabel>;

 public:
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
  template <typename DAGSRange>
  inline void AddDAGs(const DAGSRange& dags, NodeId below = {});

  template <typename DAG>
  inline void AddDAG(DAG dag, NodeId below = {});

  /**
   * Get the DAG resulting from merge
   * @{
   */
  inline MergeDAG GetResult() const;
  /** @} */

  /**
   * Access the labels of the resulting DAG's nodes.
   */
  inline const ConcurrentUnorderedMap<NodeLabel, NodeId>& GetResultNodes() const;

  inline const ConcurrentUnorderedMap<NodeId, NodeLabel>& GetResultNodeLabels() const;

  /**
   * Compute the mutations on the resulting DAG's edges and store in the result MADAG.
   */
  inline void ComputeResultEdgeMutations();

  inline bool ContainsLeafset(const LeafSet& leafset) const;

 private:
  inline MutableMergeDAG ResultDAG();

  template <typename DAGSRange>
  static void MergeCompactGenomes(size_t i, const DAGSRange& dags, NodeId below,
                                  std::vector<NodeLabelsContainer>& dags_labels,
                                  MutableMergeDAG result_dag);

  template <typename DAGSRange>
  static void ComputeLeafSets(size_t i, const DAGSRange& dags, NodeId below,
                              std::vector<NodeLabelsContainer>& dags_labels,
                              ConcurrentUnorderedSet<LeafSet>& all_leaf_sets);

  template <typename DAGSRange>
  static void MergeNodes(size_t i, const DAGSRange& dags, NodeId below,
                         std::vector<NodeLabelsContainer>& dags_labels,
                         ConcurrentUnorderedMap<NodeLabel, NodeId>& result_nodes,
                         ConcurrentUnorderedMap<NodeId, NodeLabel>& result_node_labels,
                         std::atomic<size_t>& node_id);

  template <typename DAGSRange>
  static void MergeEdges(
      size_t i, const DAGSRange& dags, NodeId below,
      std::vector<NodeLabelsContainer>& dags_labels,
      const ConcurrentUnorderedMap<NodeLabel, NodeId>& result_nodes,
      ConcurrentUnorderedMap<EdgeLabel, EdgeId>& result_edges,
      tbb::concurrent_vector<std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>>&
          added_edges);

  static inline void BuildResult(
      size_t i,
      tbb::concurrent_vector<std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>>&
          added_edges,
      std::atomic<size_t>& edge_id,
      const ConcurrentUnorderedMap<NodeLabel, NodeId>& result_nodes,
      ConcurrentUnorderedMap<EdgeLabel, EdgeId>& result_edges,
      MutableMergeDAG result_dag);

  // Every unique node leaf set, found among all input DAGs.
  ConcurrentUnorderedSet<LeafSet> all_leaf_sets_;

  // Node ids of the resulting DAG's nodes.
  ConcurrentUnorderedMap<NodeLabel, NodeId> result_nodes_;
  ConcurrentUnorderedMap<NodeId, NodeLabel> result_node_labels_;

  // Edge ids of the resulting DAG's edges.
  ConcurrentUnorderedMap<EdgeLabel, EdgeId> result_edges_;

  // Resulting DAG from merging the input DAGs.
  MergeDAGStorage<> result_dag_storage_;

  std::mutex add_dags_mtx_;
};

#include "larch/impl/merge/merge_impl.hpp"
