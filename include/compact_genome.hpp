#pragma once

#include <optional>

#include "common.hpp"

class CompactGenome {
  std::vector<std::pair<MutationPosition, char>> mutations_ = {};
  size_t hash_ = {};

 public:
  static const CompactGenome* Empty();
  CompactGenome() = default;
  CompactGenome(CompactGenome&&) = default;
  CompactGenome(const CompactGenome&) = delete;
  CompactGenome& operator=(CompactGenome&&) = default;
  CompactGenome& operator=(const CompactGenome&) = delete;

  CompactGenome(const Mutations& mutations, const CompactGenome& parent,
                std::string_view reference_sequence);

  CompactGenome(std::vector<std::pair<MutationPosition, char>>&& mutations);

  bool operator==(const CompactGenome& rhs) const noexcept;

  size_t Hash() const noexcept;

  std::optional<char> operator[](MutationPosition pos) const;

  auto begin() const -> decltype(mutations_.begin());
  auto end() const -> decltype(mutations_.end());

  bool empty() const;

  static Mutations ToEdgeMutations(std::string_view reference_sequence,
                                   const CompactGenome& parent,
                                   const CompactGenome& child);

 private:
  static size_t ComputeHash(
      const std::vector<std::pair<MutationPosition, char>>& mutations);
};

template <>
struct std::hash<CompactGenome> {
  std::size_t operator()(const CompactGenome& cg) const noexcept { return cg.Hash(); }
};

template <>
struct std::equal_to<CompactGenome> {
  std::size_t operator()(const CompactGenome& lhs,
                         const CompactGenome& rhs) const noexcept {
    return lhs == rhs;
  }
};