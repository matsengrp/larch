#pragma once

#include <cstddef>

class CompactGenome;
class LeafSet;

/*
 * NodeLabel stores the data which formally defines a node: a CompactGenome,
 * and a set of child clades, stored in a LeafSet. This allows appropriate
 * equality testing between node objects during MAD operations, like merging.
 */
class NodeLabel {
 public:
  inline NodeLabel();
  inline NodeLabel(const CompactGenome* cg, const LeafSet* ls);

  inline const CompactGenome* GetCompactGenome() const;
  inline const LeafSet* GetLeafSet() const;

  inline void SetCompactGenome(const CompactGenome* cg);
  inline void SetLeafSet(const LeafSet* ls);

  inline bool operator==(const NodeLabel& rhs) const noexcept;

  [[nodiscard]] inline size_t Hash() const noexcept;

 private:
  const CompactGenome* compact_genome_;
  const LeafSet* leaf_set_;
};

template <typename View>
class FeatureReader<NodeLabel, View> {
 public:
  const NodeLabel& GetNodeLabel();
};

template <typename View>
class FeatureWriter<NodeLabel, View> : public FeatureReader<NodeLabel, View> {
 public:
  void SetNodeLabel(NodeLabel&& node_label);
};

template <>
struct std::hash<NodeLabel> {
  inline std::size_t operator()(const NodeLabel& nl) const noexcept;
};

template <>
struct std::equal_to<NodeLabel> {
  inline bool operator()(const NodeLabel& lhs, const NodeLabel& rhs) const noexcept;
};

#include "larch/impl/merge/node_label_impl.hpp"
