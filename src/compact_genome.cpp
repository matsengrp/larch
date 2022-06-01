#include "compact_genome.hpp"

const CompactGenome* CompactGenome::Empty() {
  static CompactGenome empty = {};
  return &empty;
}

CompactGenome::CompactGenome(const Mutations& mutations, const CompactGenome& parent,
                             std::string_view reference_sequence)
    : mutations_{[&] {
        std::vector<std::pair<MutationPosition, char>> result{parent.mutations_.begin(),
                                                              parent.mutations_.end()};
        for (auto [pos, base] : mutations) {
          const bool is_valid = base != reference_sequence.at(pos.value - 1);
          auto it =
              std::lower_bound(result.begin(), result.end(), pos,
                               [](std::pair<MutationPosition, char> lhs,
                                  MutationPosition rhs) { return lhs.first < rhs; });
          if (it != result.end() and it->first == pos) {
            if (is_valid) {
              it->second = base;
            } else {
              result.erase(it);
            }
          } else {
            if (is_valid) {
              result.insert(it, {pos, base});
            }
          }
        }
        return result;
      }()},
      hash_{ComputeHash(mutations_)} {}

CompactGenome::CompactGenome(std::vector<std::pair<MutationPosition, char>>&& mutations)
    : mutations_{mutations}, hash_{ComputeHash(mutations_)} {}

bool CompactGenome::operator==(const CompactGenome& rhs) const noexcept {
  if (hash_ != rhs.hash_) return false;
  return mutations_ == rhs.mutations_;
}

bool CompactGenome::operator<(const CompactGenome& rhs) const noexcept {
  return mutations_ < rhs.mutations_;
}

size_t CompactGenome::Hash() const noexcept { return hash_; }

std::optional<char> CompactGenome::operator[](MutationPosition pos) const {
  auto it = std::lower_bound(mutations_.begin(), mutations_.end(), pos,
                             [](std::pair<MutationPosition, char> lhs,
                                MutationPosition rhs) { return lhs.first < rhs; });
  if (it != mutations_.end() and it->first == pos) {
    return it->second;
  } else {
    return std::nullopt;
  }
}

auto CompactGenome::begin() const -> decltype(mutations_.begin()) {
  return mutations_.begin();
}

auto CompactGenome::end() const -> decltype(mutations_.end()) {
  return mutations_.end();
}

bool CompactGenome::empty() const { return mutations_.empty(); }

CompactGenome CompactGenome::Copy() const {
  CompactGenome result;
  result.mutations_ = mutations_;
  result.hash_ = hash_;
  return result;
}

Mutations CompactGenome::ToEdgeMutations(std::string_view reference_sequence,
                                         const CompactGenome& parent,
                                         const CompactGenome& child) {
  Mutations result;
  for (auto [pos, child_base] : child) {
    char parent_base = reference_sequence.at(pos.value - 1);
    auto opt_parent_base = parent[pos];
    if (opt_parent_base.has_value()) {
      parent_base = opt_parent_base.value();
    }
    if (parent_base != child_base) {
      result[pos] = child_base;
    }
  }

  for (auto [pos, parent_base] : parent) {
    char child_base = reference_sequence.at(pos.value - 1);
    auto opt_child_base = child[pos];
    if (opt_child_base.has_value()) {
      child_base = opt_child_base.value();
    }
    if (child_base != parent_base) {
      result[pos] = child_base;
    }
  }
  return result;
}

size_t CompactGenome::ComputeHash(
    const std::vector<std::pair<MutationPosition, char>>& mutations) {
  size_t result = 0;
  for (auto [pos, base] : mutations) {
    result = HashCombine(result, pos.value);
    result = HashCombine(result, base);
  }
  return result;
}
