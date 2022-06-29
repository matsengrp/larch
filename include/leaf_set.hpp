#pragma once

#include "common.hpp"

class CompactGenome;
class NodeLabel;

/**
 * LeafSet provides a container which keeps track of a node's child clades.
 * A node's child clades are the sets of leaves reachable from children of
 * that node.
 * Child clades are stored in a LeafSet as a vector of vectors of leaf node
 * CompactGenomes.
 */
class LeafSet {
  std::vector<std::vector<const CompactGenome*>> clades_ = {};
  size_t hash_ = {};

 public:
  static const LeafSet* Empty();
  LeafSet() = default;
  LeafSet(LeafSet&&) = default;
  LeafSet(const LeafSet&) = delete;
  LeafSet& operator=(LeafSet&&) = default;
  LeafSet& operator=(const LeafSet&) = delete;

  LeafSet(Node node, const std::vector<NodeLabel>& labels,
          std::vector<LeafSet>& computed_leafsets);

  LeafSet(std::vector<std::vector<const CompactGenome*>>&& clades);

  bool operator==(const LeafSet& rhs) const noexcept;

  size_t Hash() const noexcept;

  auto begin() const -> decltype(clades_.begin());
  auto end() const -> decltype(clades_.end());
  bool empty() const;
  size_t size() const;

 private:
  static size_t ComputeHash(
      const std::vector<std::vector<const CompactGenome*>>& clades);
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
