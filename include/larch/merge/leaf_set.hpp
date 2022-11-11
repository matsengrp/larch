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
 public:
  inline static const LeafSet* Empty();
  LeafSet() = default;
  MOVE_ONLY(LeafSet);

  template <typename Node>
  LeafSet(Node node, const std::vector<NodeLabel>& labels,
          std::vector<LeafSet>& computed_leafsets);

  inline LeafSet(std::vector<std::vector<const CompactGenome*>>&& clades);

  inline bool operator==(const LeafSet& rhs) const noexcept;

  [[nodiscard]] inline size_t Hash() const noexcept;

  inline auto begin() const;
  inline auto end() const;
  inline bool empty() const;
  inline size_t size() const;

  [[nodiscard]] inline std::vector<const CompactGenome*> ToParentClade() const;

  inline const std::vector<std::vector<const CompactGenome*>>& GetClades() const;

 private:
  inline static size_t ComputeHash(
      const std::vector<std::vector<const CompactGenome*>>& clades) noexcept;
  std::vector<std::vector<const CompactGenome*>> clades_ = {};
  size_t hash_ = {};
};

template <typename View>
class FeatureReader<LeafSet, View> {
 public:
  const LeafSet& GetLeafSet();
};

template <typename View>
class FeatureWriter<LeafSet, View> : public FeatureReader<LeafSet, View> {
 public:
  void SetLeafSet(LeafSet&& leaf_set);
};

template <>
struct std::hash<LeafSet> {
  inline std::size_t operator()(const LeafSet& ls) const noexcept;
};

template <>
struct std::equal_to<LeafSet> {
  inline bool operator()(const LeafSet& lhs, const LeafSet& rhs) const noexcept;
};

#include "larch/impl/merge/leaf_set_impl.hpp"
