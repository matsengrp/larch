#include "larch/merge/leaf_set.hpp"

#include "larch/dag/dag.hpp"
#include "larch/merge/node_label.hpp"

const LeafSet* LeafSet::Empty() {
  static LeafSet empty = {};
  return &empty;
}

size_t LeafSet::Hash() const noexcept { return hash_; }

auto LeafSet::begin() const -> decltype(clades_.begin()) { return clades_.begin(); }

auto LeafSet::end() const -> decltype(clades_.end()) { return clades_.end(); }

bool LeafSet::empty() const { return clades_.empty(); }

size_t LeafSet::size() const { return clades_.size(); }

std::vector<const CompactGenome*> LeafSet::ToParentClade() const {
  std::vector<const CompactGenome*> result =
      ranges::to_vector(clades_ | ranges::views::join);
  result |= ranges::actions::sort | ranges::actions::unique;
  return result;
}

const std::vector<std::vector<const CompactGenome*>>& LeafSet::GetClades() const {
  return clades_;
}
