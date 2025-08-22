#pragma once

#include <map>
#include <type_traits>

#include "larch/dag/dag.hpp"
#include "larch/contiguous_map.hpp"

/**
 * @brief A wrapper for size_t, storing a 1-based index on the reference sequence.
 */
struct MutationPosition {
  size_t value = NoId;

  friend std::ostream& operator<<(std::ostream& os, MutationPosition pos) {
    os << pos.value;
    return os;
  }
};

inline bool operator==(MutationPosition lhs, MutationPosition rhs);
inline bool operator!=(MutationPosition lhs, MutationPosition rhs);
inline bool operator<(MutationPosition lhs, MutationPosition rhs);

/**
 * @brief Container for mutations associated with an edge in a mutation-annotated phylogenetic DAG.
 * 
 * EdgeMutations stores a collection of mutations that occur along a specific edge in the tree.
 * Each mutation is represented as a position on the reference sequence along with the parent
 * and child bases at that position. The class uses a ContiguousMap for efficient storage and
 * retrieval of mutations by position.
 */
class EdgeMutations {
  ContiguousMap<MutationPosition, std::pair<MutationBase, MutationBase>> mutations_;

 public:
  EdgeMutations() = default;
  MOVE_ONLY(EdgeMutations);
  template <typename T>
  EdgeMutations(T&& mutations_view,
                std::enable_if_t<not std::is_same_v<T, EdgeMutations>>* = nullptr);

  template <typename CRTP>
  [[nodiscard]] inline EdgeMutations Copy(const CRTP*) const;
  inline auto begin() const -> decltype(mutations_.begin());
  inline auto end() const -> decltype(mutations_.end());
  inline size_t size() const;
  inline bool empty() const;
  inline auto operator[](MutationPosition pos) -> decltype(mutations_[pos]);
  inline auto insert(
      std::pair<MutationPosition, std::pair<MutationBase, MutationBase>> mut);
  inline bool HasMutationAtPosition(MutationPosition pos) const;
  inline std::pair<MutationBase, MutationBase> GetMutation(MutationPosition pos) const;

  inline bool operator==(const EdgeMutations& rhs) const;
  inline bool operator!=(const EdgeMutations& rhs) const;

  inline std::string ToString() const;

 private:
  inline explicit EdgeMutations(
      ContiguousMap<MutationPosition, std::pair<MutationBase, MutationBase>>&&
          mutations);
};

template <typename CRTP, typename Tag>
struct FeatureConstView<EdgeMutations, CRTP, Tag> {
  const EdgeMutations& GetEdgeMutations() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<EdgeMutations, CRTP, Tag> {
  void SetEdgeMutations(EdgeMutations&& edge_mutations) const;
  EdgeMutations& GetMutableEdgeMutations() const;
};

#include "larch/impl/madag/edge_mutations_impl.hpp"
