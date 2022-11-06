/**
 * An EdgeLabel stores an edge's parent and child NodeLabel objects
 */

#pragma once

#include "larch/merge/node_label.hpp"

#include "larch/common.hpp"

class EdgeLabel {
 public:
  EdgeLabel() = default;
  inline EdgeLabel(NodeLabel parent, NodeLabel child);

  inline NodeLabel GetParent() const;
  inline NodeLabel GetChild() const;

  inline bool operator==(const EdgeLabel& rhs) const noexcept;

  [[nodiscard]] inline size_t Hash() const noexcept;

  [[nodiscard]] inline CladeIdx ComputeCladeIdx() const;

 private:
  NodeLabel parent_;
  NodeLabel child_;
};

template <>
struct std::hash<EdgeLabel> {
  inline std::size_t operator()(const EdgeLabel& el) const noexcept;
};

template <>
struct std::equal_to<EdgeLabel> {
  inline bool operator()(const EdgeLabel& lhs, const EdgeLabel& rhs) const noexcept;
};

#include "larch/impl/merge/edge_label_impl.hpp"
