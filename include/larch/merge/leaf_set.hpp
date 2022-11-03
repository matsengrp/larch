#pragma once

#include "larch/common.hpp"
#include "larch/madag/compact_genome.hpp"

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

  template <typename Node>
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

template <typename Node>
LeafSet::LeafSet(Node node, const std::vector<NodeLabel>& labels,
                 std::vector<LeafSet>& computed_leafsets)
    : clades_{[&] {
        std::vector<std::vector<const CompactGenome*>> clades;
        clades.reserve(node.GetCladesCount());
        for (auto clade : node.GetClades()) {
          std::vector<const CompactGenome*> clade_leafs;
          clade_leafs.reserve(clade.size());
          for (Node child : clade | Transform::GetChild()) {
            if (child.IsLeaf()) {
              clade_leafs.push_back(labels.at(child.GetId().value).GetCompactGenome());
            } else {
              for (auto& child_leafs :
                   computed_leafsets.at(child.GetId().value).clades_) {
                clade_leafs.insert(clade_leafs.end(), child_leafs.begin(),
                                   child_leafs.end());
              }
            }
          }
          clade_leafs |= ranges::actions::sort | ranges::actions::unique;
          clades.emplace_back(std::move(clade_leafs));
        }
        clades |= ranges::actions::sort;
        return clades;
      }()},
      hash_{ComputeHash(clades_)} {}

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