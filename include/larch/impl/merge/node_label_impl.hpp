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
  Assert(not id->empty());
  sample_id_ = id;
}

bool NodeLabel::operator==(const NodeLabel& rhs) const noexcept {
  return compact_genome_ == rhs.compact_genome_ && leaf_set_ == rhs.leaf_set_ &&
         sample_id_ == rhs.sample_id_;
}

size_t NodeLabel::Hash() const noexcept {
  size_t hash = HashCombine(reinterpret_cast<std::uintptr_t>(compact_genome_),
                            reinterpret_cast<std::uintptr_t>(leaf_set_));
  return HashCombine(hash, reinterpret_cast<std::uintptr_t>(sample_id_));
}

std::size_t std::hash<NodeLabel>::operator()(const NodeLabel& nl) const noexcept {
  return nl.Hash();
}

bool std::equal_to<NodeLabel>::operator()(const NodeLabel& lhs,
                                          const NodeLabel& rhs) const noexcept {
  return lhs == rhs;
}
