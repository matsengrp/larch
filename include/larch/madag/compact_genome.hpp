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
  inline static const CompactGenome* Empty();
  CompactGenome() = default;
  MOVE_ONLY(CompactGenome);
  inline CompactGenome(std::vector<std::pair<MutationPosition, char>>&& mutations);

  inline bool operator==(const CompactGenome& rhs) const noexcept;
  inline bool operator<(const CompactGenome& rhs) const noexcept;
  [[nodiscard]] inline size_t Hash() const noexcept;

  inline std::optional<char> operator[](MutationPosition pos) const;
  inline std::string ToString() const;

  inline auto begin() const -> decltype(mutations_.begin());
  inline auto end() const -> decltype(mutations_.end());

  inline bool empty() const;

  [[nodiscard]] inline CompactGenome Copy() const;

  inline void AddParentEdge(const EdgeMutations& mutations, const CompactGenome& parent,
                            std::string_view reference_sequence);

  [[nodiscard]] inline static EdgeMutations ToEdgeMutations(
      std::string_view reference_sequence, const CompactGenome& parent,
      const CompactGenome& child);

 private:
  DAG_FEATURE_FRIENDS;
  inline static size_t ComputeHash(
      const std::vector<std::pair<MutationPosition, char>>& mutations);
};

template <typename View>
class FeatureReader<CompactGenome, View> {
 public:
  const CompactGenome& GetCompactGenome();
};

template <typename View>
class FeatureWriter<CompactGenome, View> : public FeatureReader<CompactGenome, View> {
 public:
  CompactGenome& GetCompactGenome();
  void SetCompactGenome(CompactGenome&& compact_genome);
};

template <>
struct std::hash<CompactGenome> {
  inline std::size_t operator()(const CompactGenome& cg) const noexcept;
};

template <>
struct std::equal_to<CompactGenome> {
  inline bool operator()(const CompactGenome& lhs,
                         const CompactGenome& rhs) const noexcept;
};

#include "larch/impl/madag/compact_genome_impl.hpp"
