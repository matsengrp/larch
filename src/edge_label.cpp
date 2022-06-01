#include "edge_label.hpp"

#include "compact_genome.hpp"
#include "leaf_set.hpp"

bool EdgeLabel::operator==(const EdgeLabel& rhs) const noexcept {
  return parent_compact_genome == rhs.parent_compact_genome &&
         parent_leaf_set == rhs.parent_leaf_set &&
         child_compact_genome == rhs.child_compact_genome &&
         child_leaf_set == rhs.child_leaf_set;
}

size_t EdgeLabel::Hash() const noexcept {
  size_t hash = reinterpret_cast<std::uintptr_t>(parent_compact_genome);
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(parent_leaf_set));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(child_compact_genome));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(child_leaf_set));
  return hash;
}
