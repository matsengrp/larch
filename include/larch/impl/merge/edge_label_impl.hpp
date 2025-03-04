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
  auto parent_clade = child_.GetLeafSet()->ToParentClade(child_.GetSampleId());
  CladeIdx result{0};
  for (const auto& clade : *parent_.GetLeafSet()) {
    if (clade == parent_clade) {
      return result;
    }
    ++result.value;
  }
  Fail("Can't compute clade index");
  return {};
}

std::size_t std::hash<EdgeLabel>::operator()(const EdgeLabel& el) const noexcept {
  return el.Hash();
}

bool std::equal_to<EdgeLabel>::operator()(const EdgeLabel& lhs,
                                          const EdgeLabel& rhs) const noexcept {
  return lhs == rhs;
}
