#include "node_label.hpp"

#include "compact_genome.hpp"
#include "leaf_set.hpp"

NodeLabel::NodeLabel()
    : compact_genome_{CompactGenome::Empty()}, leaf_set_{LeafSet::Empty()} {}

NodeLabel::NodeLabel(const CompactGenome* cg, const LeafSet* ls)
    : compact_genome_{cg}, leaf_set_{ls} {
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

void NodeLabel::SetCompactGenome(const CompactGenome* cg) {
  compact_genome_ = cg;
  Assert(compact_genome_);
}

void NodeLabel::SetLeafSet(const LeafSet* ls) {
  leaf_set_ = ls;
  Assert(leaf_set_);
}

bool NodeLabel::operator==(const NodeLabel& rhs) const noexcept {
  return compact_genome_ == rhs.compact_genome_ && leaf_set_ == rhs.leaf_set_;
}

size_t NodeLabel::Hash() const noexcept {
  return HashCombine(reinterpret_cast<std::uintptr_t>(compact_genome_),
                     reinterpret_cast<std::uintptr_t>(leaf_set_));
}
