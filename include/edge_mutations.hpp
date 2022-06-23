#pragma once

#include <map>

#include "common.hpp"

struct MutationPosition {
  size_t value = NoId;
};

inline bool operator==(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value == rhs.value;
}

inline bool operator<(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value < rhs.value;
}

class EdgeMutations {
  std::map<MutationPosition, std::pair<char, char>> mutations_;

 public:
  EdgeMutations() = default;
  EdgeMutations(EdgeMutations&&) = default;
  EdgeMutations& operator=(EdgeMutations&&) = default;
  EdgeMutations(const EdgeMutations&) = delete;
  EdgeMutations& operator=(const EdgeMutations&) = delete;

  auto begin() const -> decltype(mutations_.begin());
  auto end() const -> decltype(mutations_.end());
  auto operator[](MutationPosition pos) -> decltype(mutations_[pos]);
  auto insert(std::pair<MutationPosition, std::pair<char, char>> mut)
      -> decltype(mutations_.insert(mut));
  bool operator==(const EdgeMutations& rhs) const;
  bool operator!=(const EdgeMutations& rhs) const;
};
