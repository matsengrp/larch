#pragma once

class CladeUnion {
 public:
  CladeUnion() = default;
  MOVE_ONLY(CladeUnion);
  //   inline explicit CladeUnion(
  //       const std::vector<const CladeUnion*>& clades);
  inline explicit CladeUnion(const std::vector<const CompactGenome*>& leafs);

  inline bool operator==(const CladeUnion& rhs) const noexcept;

  inline size_t Hash() const noexcept;

  inline auto begin() const;
  inline auto end() const;
  inline bool empty() const;
  inline size_t size() const;

 private:
  inline static size_t ComputeHash(
      const std::vector<const CompactGenome*>& leafs) noexcept;

  std::vector<const CompactGenome*> leafs_ = {};
  size_t hash_ = {};
};

template <>
struct std::hash<CladeUnion> {
  inline std::size_t operator()(const CladeUnion& cu) const noexcept;
};

template <>
struct std::equal_to<CladeUnion> {
  inline bool operator()(const CladeUnion& lhs, const CladeUnion& rhs) const noexcept;
};

#include "larch/impl/merge/clade_union_impl.hpp"