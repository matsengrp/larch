#include "edge_label.hpp"

#include "common.hpp"

EdgeLabel::EdgeLabel(NodeLabel parent, NodeLabel child)
    : parent_{parent}, child_{child} {}

NodeLabel EdgeLabel::GetParent() const { return parent_; }

NodeLabel EdgeLabel::GetChild() const { return child_; }

bool EdgeLabel::operator==(const EdgeLabel& rhs) const noexcept {
  return parent_ == rhs.parent_ && child_ == rhs.child_;
}

size_t EdgeLabel::Hash() const noexcept {
  return HashCombine(parent_.Hash(), child_.Hash());
}
