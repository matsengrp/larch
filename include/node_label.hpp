#pragma once

#include <functional>

class CompactGenome;
class LeafSet;

class NodeLabel {
 public:
  NodeLabel();
  NodeLabel(const CompactGenome* cg, const LeafSet* ls);

  bool operator==(const NodeLabel& rhs) const noexcept;

  size_t Hash() const noexcept;

  const CompactGenome* compact_genome;
  const LeafSet* leaf_set;
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
