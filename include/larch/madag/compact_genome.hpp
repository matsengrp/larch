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
#include <set>

#include "larch/dag/dag.hpp"
#include "larch/madag/edge_mutations.hpp"
#include "larch/contiguous_map.hpp"
#include "larch/contiguous_set.hpp"

class CompactGenome {
  ContiguousMap<MutationPosition, char> mutations_ = {};
  size_t hash_ = {};

 public:
  inline static const CompactGenome* Empty();
  CompactGenome() = default;
  MOVE_ONLY(CompactGenome);
  inline CompactGenome(ContiguousMap<MutationPosition, char>&& mutations);

  inline bool operator==(const CompactGenome& rhs) const noexcept;
  inline bool operator!=(const CompactGenome& rhs) const noexcept;
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

  inline void ApplyChanges(const ContiguousMap<MutationPosition, char>& changes);

  inline char GetBase(MutationPosition pos, std::string_view reference_sequence) const;

  inline ContiguousSet<MutationPosition> DifferingSites(
      const CompactGenome& other) const;

  [[nodiscard]] inline static EdgeMutations ToEdgeMutations(
      std::string_view reference_sequence, const CompactGenome& parent,
      const CompactGenome& child);

 private:
  inline static size_t ComputeHash(
      const ContiguousMap<MutationPosition, char>& mutations);
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

template <typename CRTP, typename Tag>
struct FeatureConstView<CompactGenome, CRTP, Tag> {
  const CompactGenome& GetCompactGenome() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<CompactGenome, CRTP, Tag> {
  auto& operator=(CompactGenome&& compact_genome) const;
};

#include "larch/impl/madag/compact_genome_impl.hpp"
