#include "larch/madag/compact_genome.hpp"
#include "larch/merge/leaf_set.hpp"

NodeLabel::NodeLabel()
    : compact_genome_{CompactGenome::Empty()},
      leaf_set_{LeafSet::Empty()},
      sample_id_{SampleId::Empty()} {}

NodeLabel::NodeLabel(const CompactGenome* cg, const LeafSet* ls, const SampleId* id)
    : compact_genome_{cg}, leaf_set_{ls}, sample_id_{id} {
  Assert(compact_genome_);
  Assert(leaf_set_);
  Assert(sample_id_);
}

const CompactGenome* NodeLabel::GetCompactGenome() const {
  Assert(compact_genome_);
  return compact_genome_;
}

const LeafSet* NodeLabel::GetLeafSet() const {
  Assert(leaf_set_);
  return leaf_set_;
}

const SampleId* NodeLabel::GetSampleId() const {
  Assert(sample_id_);
  return sample_id_;
}

void NodeLabel::SetCompactGenome(const CompactGenome* cg) {
  Assert(cg);
  compact_genome_ = cg;
}

void NodeLabel::SetLeafSet(const LeafSet* ls) {
  Assert(ls);
  leaf_set_ = ls;
}

void NodeLabel::SetSampleId(const SampleId* id) {
  Assert(id);
  Assert(not id->IsEmpty());
  sample_id_ = id;
}

bool NodeLabel::operator==(const NodeLabel& rhs) const noexcept {
  if (leaf_set_ != rhs.leaf_set_) {
    return false;
  }
  if (sample_id_->IsEmpty()) {
    return compact_genome_ == rhs.compact_genome_;
  } else {
    return sample_id_ == rhs.sample_id_;
  }
}

size_t NodeLabel::Hash() const noexcept {
  std::uintptr_t unique = sample_id_->IsEmpty()
                              ? reinterpret_cast<std::uintptr_t>(compact_genome_)
                              : reinterpret_cast<std::uintptr_t>(sample_id_);
  return HashCombine(unique, reinterpret_cast<std::uintptr_t>(leaf_set_));
}

std::size_t std::hash<NodeLabel>::operator()(const NodeLabel& nl) const noexcept {
  return nl.Hash();
}

bool std::equal_to<NodeLabel>::operator()(const NodeLabel& lhs,
                                          const NodeLabel& rhs) const noexcept {
  return lhs == rhs;
}
