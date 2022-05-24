#pragma once

#include <string_view>
#include <vector>
#include <map>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <shared_mutex>
#include <thread>
#include <atomic>

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_for_each.h>
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
  inline Merge(std::string_view reference_sequence,
               const std::vector<HistoryDAG>& trees,
               const std::vector<std::vector<Mutations>>& mutations,
               bool show_progress = false);

  inline void Run();

  inline HistoryDAG& GetResult();
  inline const HistoryDAG& GetResult() const;

 private:
  inline void ComputeCompactGenomes();

  inline void ComputeLeafSets();

  inline void MergeTrees();

  static inline std::vector<CompactGenome> ComputeCompactGenomes(
      const HistoryDAG& tree, const std::vector<Mutations>& edge_mutations,
      std::string_view reference_sequence);

  static inline std::vector<LeafSet> ComputeLeafSets(
      const HistoryDAG& tree, const std::vector<NodeLabel>& labels);

  std::string_view reference_sequence_;
  const std::vector<HistoryDAG>& trees_;
  const std::vector<std::vector<Mutations>>& mutations_;
  bool show_progress_;

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
