#include "edge_mutations.hpp"

auto EdgeMutations::begin() const -> decltype(mutations_.begin()) {
  return mutations_.begin();
}

auto EdgeMutations::end() const -> decltype(mutations_.end()) {
  return mutations_.end();
}

auto EdgeMutations::operator[](MutationPosition pos) -> decltype(mutations_[pos]) {
  return mutations_[pos];
}

auto EdgeMutations::insert(std::pair<MutationPosition, std::pair<char, char>> mut)
    -> decltype(mutations_.insert(mut)) {
  return mutations_.insert(mut);
}

bool EdgeMutations::operator==(const EdgeMutations& rhs) const {
  return mutations_ == rhs.mutations_;
}

bool EdgeMutations::operator!=(const EdgeMutations& rhs) const {
  return mutations_ != rhs.mutations_;
}
