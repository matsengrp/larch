#include "edge_label.hpp"

#include "dag.hpp"
#include "leaf_set.hpp"
#include "compact_genome.hpp"

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

CladeIdx EdgeLabel::ComputeCladeIdx() const {
  const auto* child_leaf_set = child_.GetLeafSet();
  Assert(child_leaf_set);
  auto parent_clade = child_leaf_set->ToParentClade();
  if (parent_clade.empty()) {
    const auto* child_compact_genome = child_.GetCompactGenome();
    Assert(child_compact_genome);
    parent_clade.push_back(child_compact_genome);
  }
  CladeIdx result{0};
  const auto* parent_leaf_set = parent_.GetLeafSet();
  Assert(parent_leaf_set);
  for (const auto& clade : *parent_leaf_set) {
    if (clade == parent_clade) {
      return result;
    }
    ++result.value;
  }
  return {};
}
