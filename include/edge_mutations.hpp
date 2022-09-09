#pragma once

#include <map>

#include "common.hpp"

/**
 * A wrapper for size_t, storing a 1-based index on the reference sequence.
 */
struct MutationPosition {
  size_t value = NoId;
};

inline bool operator==(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value == rhs.value;
}

inline bool operator<(MutationPosition lhs, MutationPosition rhs) {
  return lhs.value < rhs.value;
}

/**
 * A container for mutations along an edge, implemented as a map from
 * MutationPositions to (char, char) pairs in which the first char is the
 * site's state on the parent node, and the second is the site's state on the
 * child.
 */
class EdgeMutations {
  std::map<MutationPosition, std::pair<char, char>> mutations_;

 public:
  EdgeMutations() = default;
  EdgeMutations(EdgeMutations&&) = default;
  EdgeMutations& operator=(EdgeMutations&&) = default;
  EdgeMutations(const EdgeMutations&) = delete;
  EdgeMutations& operator=(const EdgeMutations&) = delete;
  template<typename T>
  EdgeMutations(T view):mutations_(view.begin(),view.end()){}
  [[nodiscard]] EdgeMutations Copy() const;

  auto begin() const -> decltype(mutations_.begin());
  auto end() const -> decltype(mutations_.end());
  size_t size() const;
  auto operator[](MutationPosition pos) -> decltype(mutations_[pos]);
  auto insert(std::pair<MutationPosition, std::pair<char, char>> mut)
      -> decltype(mutations_.insert(mut));
  bool operator==(const EdgeMutations& rhs) const;
  bool operator!=(const EdgeMutations& rhs) const;

 private:
  explicit EdgeMutations(
      const std::map<MutationPosition, std::pair<char, char>>& mutations);
};
