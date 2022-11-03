/**
 * A CompactGenome stores a sequence as a diff relative to a reference
 * sequence. This is implemented as a sorted vector of position, character
 * pairs. The position is a 1-based index on the reference sequence at which
 * the CompactGenome differs from that reference, and the character describes
 * that differing state.
 */
#pragma once

#include <optional>
#include <ostream>
#include <string>

#include "larch/dag/dag.hpp"
#include "larch/madag/edge_mutations.hpp"

class CompactGenome {
  std::vector<std::pair<MutationPosition, char>> mutations_ = {};
  size_t hash_ = {};

 public:
  static const CompactGenome* Empty();
  CompactGenome() = default;
  MOVE_ONLY(CompactGenome);
  CompactGenome(std::vector<std::pair<MutationPosition, char>>&& mutations);

  bool operator==(const CompactGenome& rhs) const noexcept;
  bool operator<(const CompactGenome& rhs) const noexcept;
  [[nodiscard]] size_t Hash() const noexcept;

  std::optional<char> operator[](MutationPosition pos) const;
  std::string ToString() const;

  auto begin() const -> decltype(mutations_.begin());
  auto end() const -> decltype(mutations_.end());

  bool empty() const;

  [[nodiscard]] CompactGenome Copy() const;

  void AddParentEdge(const EdgeMutations& mutations, const CompactGenome& parent,
                     std::string_view reference_sequence);

  [[nodiscard]] static EdgeMutations ToEdgeMutations(
      std::string_view reference_sequence, const CompactGenome& parent,
      const CompactGenome& child);

 private:
  DAG_FEATURE_FRIENDS;
  static size_t ComputeHash(
      const std::vector<std::pair<MutationPosition, char>>& mutations);
};

template <typename View>
class FeatureReader<CompactGenome, View> {
 public:
  const CompactGenome& GetCompactGenome() const;
};

template <typename View>
class FeatureWriter<CompactGenome, View> : public FeatureReader<CompactGenome, View> {
 public:
   void SetCompactGenome(CompactGenome&& compact_genome);
};

template <>
struct std::hash<CompactGenome> {
  std::size_t operator()(const CompactGenome& cg) const noexcept;
};

template <>
struct std::equal_to<CompactGenome> {
  std::size_t operator()(const CompactGenome& lhs,
                         const CompactGenome& rhs) const noexcept;
};

#include "larch/impl/madag/compact_genome_impl.hpp"