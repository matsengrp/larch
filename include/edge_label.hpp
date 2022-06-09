#pragma once

#include "node_label.hpp"

class EdgeLabel {
 public:
  EdgeLabel() = default;
  EdgeLabel(NodeLabel parent, NodeLabel child);

  NodeLabel GetParent() const;
  NodeLabel GetChild() const;

  bool operator==(const EdgeLabel& rhs) const noexcept;

  size_t Hash() const noexcept;

 private:
  NodeLabel parent_;
  NodeLabel child_;
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
