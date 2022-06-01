#include "node_label.hpp"

#include "compact_genome.hpp"
#include "leaf_set.hpp"

NodeLabel::NodeLabel()
    : compact_genome{CompactGenome::Empty()}, leaf_set{LeafSet::Empty()} {}

NodeLabel::NodeLabel(const CompactGenome* cg, const LeafSet* ls)
    : compact_genome{cg}, leaf_set{ls} {}

bool NodeLabel::operator==(const NodeLabel& rhs) const noexcept {
  return compact_genome == rhs.compact_genome && leaf_set == rhs.leaf_set;
}

size_t NodeLabel::Hash() const noexcept {
  return HashCombine(reinterpret_cast<std::uintptr_t>(compact_genome),
                     reinterpret_cast<std::uintptr_t>(leaf_set));
}
