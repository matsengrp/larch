#pragma once

#include "larch/common.hpp"
#include "larch/madag/sample_id.hpp"

class NodeLabel;

/**
 * LeafSet provides a container which keeps track of a node's child clades.
 * A node's child clades are the sets of leaves reachable from children of
 * that node.
 * Child clades are stored in a LeafSet as a vector of vectors of leaf node
 * CompactGenomes.
 */
class LeafSet {
  using UniqueData = const SampleId*;
  std::vector<std::vector<UniqueData>> clades_ = {};
  size_t hash_ = {};

 public:
  inline static const LeafSet* Empty();
  LeafSet() = default;
  MOVE_ONLY(LeafSet);

  template <typename Node>
  LeafSet(Node node, const std::vector<NodeLabel>& labels,
          std::vector<LeafSet>& computed_leafsets);

  inline LeafSet(std::vector<std::vector<UniqueData>>&& clades);

  inline bool operator==(const LeafSet& rhs) const noexcept;

  [[nodiscard]] inline size_t Hash() const noexcept;

  inline auto begin() const -> decltype(clades_.begin());
  inline auto end() const -> decltype(clades_.end());
  inline bool empty() const;
  inline size_t size() const;

  [[nodiscard]] inline std::vector<UniqueData> ToParentClade() const;

  [[nodiscard]] inline size_t ParentCladeSize() const;

  inline const std::vector<std::vector<UniqueData>>& GetClades() const;

  inline std::string ToString() const;

  template <typename DAGType>
  inline static std::vector<LeafSet> ComputeLeafSets(
      DAGType dag, const std::vector<NodeLabel>& labels);

 private:
  inline static size_t ComputeHash(
      const std::vector<std::vector<UniqueData>>& clades) noexcept;
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
