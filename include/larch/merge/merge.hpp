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

using MergeDAGStorage =
    ExtendDAGStorage<DefaultDAGStorage,
                     Extend::Nodes<Deduplicate<CompactGenome>, SampleId>,
                     Extend::Edges<EdgeMutations>, Extend::DAG<ReferenceSequence>>;

using MergeDAG = DAGView<const MergeDAGStorage>;
using MutableMergeDAG = DAGView<MergeDAGStorage>;

template <typename DAG>
class Fragment {
 public:
  template <typename Id, typename Feature>
  static const bool contains_element_feature =
      DAG::template contains_element_feature<Id, Feature>;

  Fragment(DAG dag, std::vector<NodeId>&& nodes, std::vector<EdgeId> edges)
      : dag_{dag}, nodes_{nodes}, edges_{edges} {}

  void AssertUA() const { dag_.AssertUA(); }

  size_t GetNodesCount() const { return nodes_.size(); }

  size_t GetEdgesCount() const { return edges_.size(); }

  auto Get(NodeId id) const { return dag_.Get(id); }

  auto Get(EdgeId id) const { return dag_.Get(id); }

  auto GetNodes() const {
    return ranges::views::all(nodes_) | Transform::ToNodes(dag_);
  }

  auto GetEdges() const {
    return ranges::views::all(edges_) | Transform::ToEdges(dag_);
  }

  auto GetRoot() const { return dag_.GetRoot(); }

 private:
  DAG dag_;
  const std::vector<NodeId> nodes_;
  const std::vector<EdgeId> edges_;
};

class Merge {
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

  /**
   * Get the DAG resulting from merge
   * @{
   */
  inline MergeDAG GetResult() const;
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
  inline MutableMergeDAG ResultDAG();

  template <typename DAGSRange>
  void ComputeLeafSets(const DAGSRange& dags,
                       std::vector<std::vector<NodeLabel>>& dags_labels, NodeId below);

  template <typename DAGSRange>
  void MergeDAGs(const DAGSRange& dags,
                 const std::vector<std::vector<NodeLabel>>& dags_labels, NodeId below);

  // Every unique node leaf set, found among all input DAGs.
  ConcurrentUnorderedSet<LeafSet> all_leaf_sets_;

  // Node ids of the resulting DAG's nodes.
  std::unordered_map<NodeLabel, NodeId> result_nodes_;
  std::vector<NodeLabel> result_node_labels_;

  // Edge ids of the resulting DAG's edges.
  ConcurrentUnorderedMap<EdgeLabel, EdgeId> result_edges_;

  // Resulting DAG from merging the input DAGs.
  MergeDAGStorage result_dag_storage_;
};

#include "larch/impl/merge/merge_impl.hpp"
