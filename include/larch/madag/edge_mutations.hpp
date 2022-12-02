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

  std::string ToString() const {
    std::string result = "<";
    for (auto [pos, muts] : mutations_) {
      result += muts.first;
      result += std::to_string(pos.value);
      result += muts.second;
      result += ",";
    }
    result += ">";
    return result;
  }

 private:
  inline explicit EdgeMutations(
      std::map<MutationPosition, std::pair<char, char>>&& mutations);
};

template <typename CRTP, typename Tag>
struct FeatureConstView<EdgeMutations, CRTP, Tag> {
  const EdgeMutations& GetEdgeMutations() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<EdgeMutations, CRTP, Tag> {
  void SetEdgeMutations(EdgeMutations&& edge_mutations) const;
};

#include "larch/impl/madag/edge_mutations_impl.hpp"
