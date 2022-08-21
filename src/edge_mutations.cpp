#include "edge_mutations.hpp"

EdgeMutations::EdgeMutations(
    const std::map<MutationPosition, std::pair<char, char>>& mutations)
    : mutations_{mutations} {}

EdgeMutations EdgeMutations::Copy() const { return EdgeMutations{mutations_}; }

auto EdgeMutations::begin() const -> decltype(mutations_.begin()) {
  return mutations_.begin();
}

auto EdgeMutations::end() const -> decltype(mutations_.end()) {
  return mutations_.end();
}

size_t EdgeMutations::size() const { return mutations_.size(); }

auto EdgeMutations::operator[](MutationPosition pos) -> decltype(mutations_[pos]) {
  Assert(pos.value != NoId);
  return mutations_[pos];
}

auto EdgeMutations::insert(std::pair<MutationPosition, std::pair<char, char>> mut)
    -> decltype(mutations_.insert(mut)) {
  Assert(mut.first.value != NoId);
  return mutations_.insert(mut);
}

bool EdgeMutations::operator==(const EdgeMutations& rhs) const {
  return mutations_ == rhs.mutations_;
}

bool EdgeMutations::operator!=(const EdgeMutations& rhs) const {
  return mutations_ != rhs.mutations_;
}
