#pragma once

#include <cstddef>

class CompactGenome;
class LeafSet;

class EdgeLabel {
 public:
  bool operator==(const EdgeLabel& rhs) const noexcept;

  size_t Hash() const noexcept;

  const CompactGenome* parent_compact_genome = nullptr;
  const LeafSet* parent_leaf_set = nullptr;
  const CompactGenome* child_compact_genome = nullptr;
  const LeafSet* child_leaf_set = nullptr;
};

namespace std {
template <typename>
struct hash;
template <typename>
struct equal_to;
}  // namespace std

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
