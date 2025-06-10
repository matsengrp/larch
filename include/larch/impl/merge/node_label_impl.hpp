#include "larch/madag/compact_genome.hpp"

NodeLabel::NodeLabel()
    : compact_genome_{CompactGenome::GetEmpty()},
      leaf_set_{LeafSet::GetEmpty()},
      sample_id_{SampleId::GetEmpty()} {}

NodeLabel::NodeLabel(const CompactGenome* cg, const LeafSet* ls, UniqueData id)
    : compact_genome_{cg}, leaf_set_{ls}, sample_id_{id} {
  Assert(compact_genome_);
  Assert(leaf_set_);
}

const CompactGenome* NodeLabel::GetCompactGenome() const {
  Assert(compact_genome_);
  return compact_genome_;
}

const LeafSet* NodeLabel::GetLeafSet() const {
  Assert(leaf_set_);
  return leaf_set_;
}

UniqueData NodeLabel::GetSampleId() const { return sample_id_; }

void NodeLabel::SetCompactGenome(const CompactGenome* cg) {
  Assert(cg);
  Assert(not cg->empty());
  Assert(compact_genome_->empty());
  Assert(sample_id_.empty());
  compact_genome_ = cg;
}

void NodeLabel::SetLeafSet(const LeafSet* ls) {
  Assert(ls);
  Assert(not ls->empty());
  Assert(leaf_set_->empty());
  leaf_set_ = ls;
}

void NodeLabel::SetSampleId(UniqueData id) {
  Assert(not id.empty());
  Assert(sample_id_.empty());
  sample_id_ = id;
}

bool NodeLabel::operator==(const NodeLabel& rhs) const noexcept {
  if (leaf_set_ != rhs.leaf_set_) {
    return false;
  }
  if (sample_id_.empty()) {
    return compact_genome_ == rhs.compact_genome_;
  } else {
    return std::equal_to<SampleId>{}(sample_id_, rhs.sample_id_);
  }
}

size_t NodeLabel::Hash() const noexcept {
  std::uintptr_t unique = sample_id_.empty()
                              ? reinterpret_cast<std::uintptr_t>(compact_genome_)
                              : static_cast<std::uintptr_t>(sample_id_.Hash());
  return HashCombine(unique, reinterpret_cast<std::uintptr_t>(leaf_set_));
}

bool NodeLabel::empty() const {
  Assert(compact_genome_);
  Assert(leaf_set_);
  if (compact_genome_->empty() and sample_id_.empty() and leaf_set_->empty()) {
    return true;
  }
  return false;
}

std::string NodeLabel::ToString() const {
  std::string result = "[";
  result += compact_genome_->ToString();
  result += "(";
  result += leaf_set_->ToString();
  result += ")";
  return result;
}

std::size_t std::hash<NodeLabel>::operator()(const NodeLabel& nl) const noexcept {
  return nl.Hash();
}

bool std::equal_to<NodeLabel>::operator()(const NodeLabel& lhs,
                                          const NodeLabel& rhs) const noexcept {
  return lhs == rhs;
}
