#pragma once

#include <string_view>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <shared_mutex>
#include <thread>
#include <atomic>
#include <numeric>

#include "larch/madag/mutation_annotated_dag.hpp"
#include "larch/merge/leaf_set.hpp"
#include "larch/merge/node_label.hpp"
#include "larch/merge/edge_label.hpp"

using MergeDAGStorage =
    ExtendDAGStorage<DefaultDAGStorage,
                     Extend::Nodes<Deduplicate<CompactGenome>, SampleId>,
                     Extend::Edges<EdgeMutations>, Extend::DAG<ReferenceSequence>>;

using MergeDAG = DAGView<const MergeDAGStorage>;
using MutableMergeDAG = DAGView<MergeDAGStorage>;

class Merge {
  using AddedEdge = std::tuple<EdgeLabel, EdgeId, NodeId, NodeId, CladeIdx>;

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
  inline void AddDAGs(DAGSRange&& dags, NodeId below = {});

  template <typename DAG>
  inline void AddDAG(DAG&& dag, NodeId below = {});

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

  template <typename DAG>
  void MergeCompactGenomes(DAG dag, NodeId below, std::vector<NodeLabel>& labels);

  template <typename DAG>
  void ComputeLeafSets(DAG dag, NodeId below, std::vector<NodeLabel>& labels);

  template <typename DAG>
  void MergeNodes(size_t worker, DAG dag, NodeId below,
                  const std::vector<NodeLabel>& labels, std::atomic<size_t>& node_id,
                  Reduction<NodeId>& added_nodes);

  template <typename DAG>
  void MergeEdges(size_t worker, DAG dag, NodeId below,
                  const std::vector<NodeLabel>& labels,
                  Reduction<AddedEdge>& added_edges);

  inline void BuildResult(AddedEdge& added_edge, std::atomic<size_t>& edge_id);

  // Every unique node leaf set, found among all input DAGs.
  ConcurrentUnorderedSet<LeafSet> all_leaf_sets_;

  // Node ids of the resulting DAG's nodes.
  ConcurrentUnorderedMap<NodeLabel, NodeId> result_nodes_;

  ConcurrentUnorderedMap<NodeId, NodeLabel> result_node_labels_;

  // Edge ids of the resulting DAG's edges.
  ConcurrentUnorderedMap<EdgeLabel, EdgeId> result_edges_;

  // Resulting DAG from merging the input DAGs.
  MergeDAGStorage result_dag_storage_;
};

#include "larch/impl/merge/merge_impl.hpp"
