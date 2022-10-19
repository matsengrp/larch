#pragma once

#include "common.hpp"
#include "compact_genome.hpp"

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

  inline LeafSet(std::vector<std::vector<const CompactGenome*>>&& clades);

  inline bool operator==(const LeafSet& rhs) const noexcept;

  [[nodiscard]] size_t Hash() const noexcept;

  auto begin() const -> decltype(clades_.begin());
  auto end() const -> decltype(clades_.end());
  bool empty() const;
  size_t size() const;

  [[nodiscard]] std::vector<const CompactGenome*> ToParentClade() const;

  const std::vector<std::vector<const CompactGenome*>>& GetClades() const;

 private:
  inline static size_t ComputeHash(
      const std::vector<std::vector<const CompactGenome*>>& clades) noexcept;
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

bool LeafSet::operator==(const LeafSet& rhs) const noexcept {
  return clades_ == rhs.clades_;
}

LeafSet::LeafSet(std::vector<std::vector<const CompactGenome*>>&& clades)
    : clades_{std::forward<std::vector<std::vector<const CompactGenome*>>>(clades)},
      hash_{ComputeHash(clades_)} {}

size_t LeafSet::ComputeHash(
    const std::vector<std::vector<const CompactGenome*>>& clades) noexcept {
  size_t hash = 0;
  for (auto& clade : clades) {
    for (auto leaf : clade) {
      hash = HashCombine(hash, leaf->Hash());
    }
  }
  return hash;
}