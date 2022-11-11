#pragma once

#include <map>
#include <type_traits>

#include "larch/dag/dag.hpp"

/**
 * A wrapper for size_t, storing a 1-based index on the reference sequence.
 */
struct MutationPosition {
  size_t value = NoId;
};

inline bool operator==(MutationPosition lhs, MutationPosition rhs);
inline bool operator<(MutationPosition lhs, MutationPosition rhs);

class EdgeMutations {
  std::map<MutationPosition, std::pair<char, char>> mutations_;

 public:
  EdgeMutations() = default;
  MOVE_ONLY(EdgeMutations);
  template <typename T>
  EdgeMutations(T&& mutations_view,
                std::enable_if_t<not std::is_same_v<T, EdgeMutations>>* = nullptr);

  [[nodiscard]] inline EdgeMutations Copy() const;
  inline auto begin() const -> decltype(mutations_.begin());
  inline auto end() const -> decltype(mutations_.end());
  inline size_t size() const;
  inline auto operator[](MutationPosition pos) -> decltype(mutations_[pos]);
  inline auto insert(std::pair<MutationPosition, std::pair<char, char>> mut)
      -> decltype(mutations_.insert(mut));
  inline bool operator==(const EdgeMutations& rhs) const;
  inline bool operator!=(const EdgeMutations& rhs) const;

 private:
  DAG_FEATURE_FRIENDS;
  inline explicit EdgeMutations(
      std::map<MutationPosition, std::pair<char, char>>&& mutations);
};

template <typename View>
class FeatureReader<EdgeMutations, View> {
 public:
  const EdgeMutations& GetEdgeMutations();
};

template <typename View>
class FeatureWriter<EdgeMutations, View> : public FeatureReader<EdgeMutations, View> {
 public:
  void SetEdgeMutations(EdgeMutations&& edge_mutations);
};

#include "larch/impl/madag/edge_mutations_impl.hpp"
