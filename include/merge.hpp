#pragma once

#include <string_view>
#include <vector>
#include <map>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <execution>
#include <shared_mutex>
#include <thread>
#include <atomic>

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>
#include <range/v3/view/enumerate.hpp>

#include "history_dag.hpp"

class CompactGenome {
 public:
  static inline const CompactGenome* Empty();
  CompactGenome() = default;
  CompactGenome(CompactGenome&&) = default;
  CompactGenome(const CompactGenome&) = delete;
  CompactGenome& operator=(CompactGenome&&) = default;
  CompactGenome& operator=(const CompactGenome&) = delete;

  inline CompactGenome(const Mutations& mutations, const CompactGenome& parent,
                       std::string_view reference_sequence);

  inline bool operator==(const CompactGenome& rhs) const;

  inline bool operator<(const CompactGenome& rhs) const;

  inline size_t Hash() const;

 private:
  std::map<MutationPosition, char> mutations_ = {};
  size_t hash_ = {};
};

class NodeLabel;

class LeafSet {
 public:
  static inline const LeafSet* Empty();
  LeafSet() = default;
  LeafSet(LeafSet&&) = default;
  LeafSet(const LeafSet&) = delete;
  LeafSet& operator=(LeafSet&&) = default;
  LeafSet& operator=(const LeafSet&) = delete;

  inline LeafSet(Node node, const std::vector<NodeLabel>& labels,
                 std::vector<LeafSet>& computed_leafsets);

  inline bool operator==(const LeafSet& rhs) const;

  inline bool operator<(const LeafSet& rhs) const;

  inline size_t Hash() const;

 private:
  std::set<std::set<const CompactGenome*>> clades_ = {};
  size_t hash_ = {};
};

class NodeLabel {
 public:
  inline NodeLabel();
  inline NodeLabel(const CompactGenome* cg, const LeafSet* ls);

  inline bool operator==(const NodeLabel& rhs) const;

  inline size_t Hash() const;

  const CompactGenome* compact_genome;
  const LeafSet* leaf_set;
};

class EdgeLabel {
 public:
  inline bool operator==(const EdgeLabel& rhs) const;

  inline size_t Hash() const;

  const CompactGenome* parent_compact_genome = nullptr;
  const LeafSet* parent_leaf_set = nullptr;
  const CompactGenome* child_compact_genome = nullptr;
  const LeafSet* child_leaf_set = nullptr;
};

template <>
struct std::hash<CompactGenome> {
  std::size_t operator()(const CompactGenome& cg) const noexcept { return cg.Hash(); }
};

template <>
struct std::equal_to<CompactGenome> {
  std::size_t operator()(const CompactGenome& lhs,
                         const CompactGenome& rhs) const noexcept {
    return lhs == rhs;
  }
};

template <>
struct std::hash<LeafSet> {
  std::size_t operator()(const LeafSet& ls) const noexcept { return ls.Hash(); }
};

template <>
struct std::equal_to<LeafSet> {
  std::size_t operator()(const LeafSet& lhs, const LeafSet& rhs) const noexcept {
    return lhs == rhs;
  }
};

template <>
struct std::hash<NodeLabel> {
  std::size_t operator()(const NodeLabel& nl) const noexcept { return nl.Hash(); }
};

template <>
struct std::equal_to<NodeLabel> {
  std::size_t operator()(const NodeLabel& lhs, const NodeLabel& rhs) const noexcept {
    return lhs == rhs;
  }
};

template <>
struct std::hash<EdgeLabel> {
  std::size_t operator()(const EdgeLabel& el) const noexcept { return el.Hash(); }
};

template <>
struct std::equal_to<EdgeLabel> {
  std::size_t operator()(const EdgeLabel& lhs, const EdgeLabel& rhs) const noexcept {
    return lhs == rhs;
  }
};

class Merge {
 public:
  Merge(std::string_view reference_sequence, const std::vector<HistoryDAG>& trees,
        const std::vector<std::vector<Mutations>>& mutations)
      : reference_sequence_{reference_sequence}, trees_{trees}, mutations_{mutations} {
    tree_labels_.resize(trees_.size());
  }

  void Run() {
    ComputeCompactGenomes();
    ComputeLeafSets();
    MergeTrees();
    std::cout << "Nodes: " << result_nodes_.size() << "\n";
    std::cout << "Edges: " << result_edges_.size() << "\n";
    result_.BuildConnections();
  }

  HistoryDAG& GetResult() { return result_; }
  const HistoryDAG& GetResult() const { return result_; }

 private:
  void ComputeCompactGenomes() {
    std::vector<size_t> tree_idxs;
    tree_idxs.resize(trees_.size());
    std::iota(tree_idxs.begin(), tree_idxs.end(), 0);

    std::cout << "Computing compact genomes " << std::flush;
    std::for_each(
        std::execution::par, tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
          const HistoryDAG& tree = trees_.at(tree_idx);
          const std::vector<Mutations>& edge_mutations = mutations_.at(tree_idx);
          std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
          labels.resize(tree.GetNodes().size());
          std::vector<CompactGenome> computed_cgs =
              ComputeCompactGenomes(tree, edge_mutations, reference_sequence_);
          for (size_t node_idx = 0; node_idx < tree.GetNodes().size(); ++node_idx) {
            auto cg_iter =
                all_compact_genomes_.insert(std::move(computed_cgs.at(node_idx)));
            labels.at(node_idx).compact_genome = std::addressof(*cg_iter.first);
          }
          std::cout << "." << std::flush;
        });
    std::cout << " done.\n";
  }

  void ComputeLeafSets() {
    std::vector<size_t> tree_idxs;
    tree_idxs.resize(trees_.size());
    std::iota(tree_idxs.begin(), tree_idxs.end(), 0);

    std::cout << "Computing leaf sets " << std::flush;

    std::for_each(
        std::execution::par, tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
          const HistoryDAG& tree = trees_.at(tree_idx);
          std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
          std::vector<LeafSet> computed_ls = ComputeLeafSets(tree, labels);
          for (size_t node_idx = 0; node_idx < tree.GetNodes().size(); ++node_idx) {
            auto ls_iter = all_leaf_sets_.insert(std::move(computed_ls.at(node_idx)));
            labels.at(node_idx).leaf_set = std::addressof(*ls_iter.first);
          }
          std::cout << "." << std::flush;
        });
    std::cout << " done.\n";
  }

  void MergeTrees() {
    std::vector<size_t> tree_idxs;
    tree_idxs.resize(trees_.size());
    std::iota(tree_idxs.begin(), tree_idxs.end(), 0);

    std::atomic<size_t> node_id{0};
    std::mutex mtx;
    std::for_each(std::execution::par, tree_idxs.begin(), tree_idxs.end(),
                  [&](size_t tree_idx) {
                    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
                    for (auto label : labels) {
                      std::lock_guard<std::mutex> lock{mtx};
                      auto i = result_nodes_.find(label);
                      if (i == result_nodes_.end()) {
                        result_nodes_.insert({label, {node_id.fetch_add(1)}});
                      }
                    }
                  });
    std::for_each(
        std::execution::par, tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
          const HistoryDAG& tree = trees_.at(tree_idx);
          const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
          for (Edge edge : tree.GetEdges()) {
            auto& parent_label = labels.at(edge.GetParent().GetId().value);
            auto& child_label = labels.at(edge.GetChild().GetId().value);
            result_edges_.insert({parent_label.compact_genome, parent_label.leaf_set,
                                  child_label.compact_genome, child_label.leaf_set});
          }
        });
    result_.InitializeComponents(result_nodes_.size(), result_edges_.size());
    size_t edge_id = 0;
    for (auto& edge : result_edges_) {
      auto parent = result_nodes_.find(
          NodeLabel{edge.parent_compact_genome, edge.parent_leaf_set});
      auto child =
          result_nodes_.find(NodeLabel{edge.child_compact_genome, edge.child_leaf_set});
      Assert(parent != result_nodes_.end());
      Assert(child != result_nodes_.end());
      Assert(parent->second.value != NoId);
      Assert(child->second.value != NoId);
      result_.AddEdge({edge_id++}, parent->second, child->second, {0});
    }
  }

  static std::vector<CompactGenome> ComputeCompactGenomes(
      const HistoryDAG& tree, const std::vector<Mutations>& edge_mutations,
      std::string_view reference_sequence) {
    std::vector<CompactGenome> result;
    result.resize(tree.GetNodes().size());
    for (auto iter : tree.TraversePreOrder()) {
      const Mutations& mutations = edge_mutations.at(iter.GetEdge().GetId().value);
      const CompactGenome& parent = result.at(iter.GetEdge().GetParent().GetId().value);
      CompactGenome& compact_genome = result.at(iter.GetNode().GetId().value);
      compact_genome = CompactGenome{mutations, parent, reference_sequence};
    }
    return result;
  }

  static std::vector<LeafSet> ComputeLeafSets(const HistoryDAG& tree,
                                              const std::vector<NodeLabel>& labels) {
    std::vector<LeafSet> result;
    result.resize(tree.GetNodes().size());
    for (Node node : tree.TraversePostOrder()) {
      result.at(node.GetId().value) = LeafSet{node, labels, result};
    }
    return result;
  }

  std::string_view reference_sequence_;
  const std::vector<HistoryDAG>& trees_;
  const std::vector<std::vector<Mutations>>& mutations_;

  tbb::concurrent_unordered_set<CompactGenome, std::hash<CompactGenome>,
                                std::equal_to<CompactGenome>>
      all_compact_genomes_;
  tbb::concurrent_unordered_set<LeafSet, std::hash<LeafSet>, std::equal_to<LeafSet>>
      all_leaf_sets_;
  std::vector<std::vector<NodeLabel>> tree_labels_;

  tbb::concurrent_unordered_map<NodeLabel, NodeId, std::hash<NodeLabel>,
                                std::equal_to<NodeLabel>>
      result_nodes_;
  tbb::concurrent_unordered_set<EdgeLabel, std::hash<EdgeLabel>,
                                std::equal_to<EdgeLabel>>
      result_edges_;
  HistoryDAG result_;
};

#include "impl/merge_impl.hpp"
