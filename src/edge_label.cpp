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
  auto parent_clade = child_.GetLeafSet()->ToParentClade();
  if (parent_clade.empty()) {
    parent_clade.push_back(child_.GetCompactGenome());
  }
  CladeIdx result{0};
  for (const auto& clade : *parent_.GetLeafSet()) {
    if (clade == parent_clade) {
      return result;
    }
    ++result.value;
  }
  return {};
}
