#pragma once

#include <vector>
#include <functional>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <cassert>
#include <string_view>
#include <algorithm>

#include "history_dag.hpp"

using CompactGenome = std::map<size_t, char>;

struct CompactGenomePointer {
  const CompactGenome* value = nullptr;
  mutable size_t hash = 0;
};

inline bool operator==(const CompactGenomePointer& lhs,
                       const CompactGenomePointer& rhs) {
  return *lhs.value == *rhs.value;
}

template <>
struct std::hash<CompactGenomePointer> {
  std::size_t operator()(const CompactGenomePointer& cg) const noexcept {
    if (cg.hash != 0) {
      return cg.hash;
    }
    size_t hash = 0;
    for (auto [pos, mut] : *cg.value) {
      hash = HashCombine(hash, pos);
      hash = HashCombine(hash, mut);
    }
    cg.hash = hash;
    return hash;
  }
};

using LeafSet = std::vector<std::unordered_set<CompactGenomePointer>>;

struct LeafSetPointer {
  const LeafSet* value = nullptr;
  mutable size_t hash = 0;
};

inline bool operator==(const LeafSetPointer& lhs, const LeafSetPointer& rhs) {
  return *lhs.value == *rhs.value;
}

template <>
struct std::hash<LeafSetPointer> {
  std::size_t operator()(const LeafSetPointer& ls) const noexcept {
    if (ls.hash != 0) {
      return ls.hash;
    }
    size_t hash = 0;
    for (auto& clade : *ls.value) {
      std::vector<CompactGenomePointer> leafs;
      for (auto& leaf : clade) {
        leafs.push_back(leaf);
      }
      std::sort(leafs.begin(), leafs.end(),
                [](const CompactGenomePointer& lhs, const CompactGenomePointer& rhs) {
                  return *lhs.value < *rhs.value;
                });
      for (auto& leaf : leafs) {
        hash = HashCombine(hash, std::hash<CompactGenomePointer>{}(leaf));
      }
    }
    ls.hash = hash;
    return hash;
  }
};

using NodeLabel = std::pair<CompactGenomePointer, LeafSetPointer>;
using EdgeLabel = std::tuple<NodeLabel, NodeLabel, CladeIdx>;

inline bool operator==(const NodeLabel& lhs, const NodeLabel& rhs) {
  return *lhs.first.value == *rhs.first.value && *lhs.second.value == *rhs.second.value;
}

template <>
struct std::hash<NodeLabel> {
  std::size_t operator()(const NodeLabel& label) const noexcept {
    size_t hash = 0;
    hash = HashCombine(hash, std::hash<CompactGenomePointer>{}(label.first));
    hash = HashCombine(hash, std::hash<LeafSetPointer>{}(label.second));
    return hash;
  }
};

inline bool operator==(const EdgeLabel& lhs, const EdgeLabel& rhs) {
  return std::get<0>(lhs) == std::get<0>(rhs) && std::get<1>(lhs) == std::get<1>(rhs);
}

template <>
struct std::hash<EdgeLabel> {
  std::size_t operator()(const EdgeLabel& label) const noexcept {
    size_t hash = 0;
    hash = HashCombine(hash, std::hash<NodeLabel>{}(std::get<0>(label)));
    hash = HashCombine(hash, std::hash<NodeLabel>{}(std::get<1>(label)));
    return hash;
  }
};

struct TreeLabels {
  std::vector<CompactGenome> compact_genomes;
  std::vector<LeafSet> leaf_sets;
  NodeLabel GetLabel(NodeId node) const {
    return {{&compact_genomes.at(node.value)}, {&leaf_sets.at(node.value)}};
  }
};

class Merge {
 public:
  inline Merge(const std::string_view& refseq,
               std::vector<std::reference_wrapper<const HistoryDAG>>&& trees,
               const std::vector<TreeLabels>& labels);

  inline void Run();

  inline HistoryDAG& GetResult();
  inline const HistoryDAG& GetResult() const;

  inline NodeLabel GetNodeLabel(size_t tree_idx, NodeId node_id);
  inline NodeId GetResultNode(NodeLabel label);
  inline void MakeResultEdge(size_t tree_idx, EdgeId edge_id);

 private:
  const std::string refseq_;
  std::vector<std::reference_wrapper<const HistoryDAG>> trees_;
  const std::vector<TreeLabels>& labels_;
  std::unordered_map<NodeLabel, NodeId> result_nodes_;
  std::unordered_set<EdgeLabel> result_edges_;
  HistoryDAG result_;
};

inline TreeLabels GetLabels(const HistoryDAG& tree, const std::string_view& refseq,
                            const std::vector<CompactGenome>& mutations);

#include "impl/merge_impl.hpp"
