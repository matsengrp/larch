/**
 * MutationBase represents a DNA base value in the MADAG. Does not support ambiguity or
 * gap characters. Bases are encoded as a two-bit array. Bases are ordered such that
 * complementary bases are represented by inverting the bits.
 */

#pragma once

#include <string_view>
#include <vector>
#include <map>

struct MutationBase {
  static constexpr size_t BitCount = 4;
  using BitArray = std::array<bool, BitCount>;

  MutationBase() = default;
  inline MutationBase(const BitArray m_value_in);
  inline MutationBase(const char m_char_in);

  inline bool IsCompatible(const MutationBase &rhs) const;
  inline bool IsAmbiguous() const;

  inline MutationBase GetComplementaryBase() const;
  inline MutationBase GetFirstBase() const;
  inline MutationBase GetCommonBases(const MutationBase &rhs) const;
  inline MutationBase GetFirstCommonBase(const MutationBase &rhs) const;

  inline char ToChar() const;
  static inline std::string ToString(std::vector<MutationBase> m_in);
  friend std::ostream &operator<<(std::ostream &os, const MutationBase m_in);
  friend std::ostream &operator<<(std::ostream &os,
                                  const std::vector<MutationBase> m_in);
  friend std::string &operator+=(std::string &str, const MutationBase m_in);

  inline bool operator==(const MutationBase &rhs) const;
  inline bool operator!=(const MutationBase &rhs) const;
  inline bool operator<(const MutationBase &rhs) const;
  inline bool operator==(const BitArray &rhs) const;
  inline bool operator!=(const BitArray &rhs) const;
  inline bool operator<(const BitArray &rhs) const;
  inline bool operator==(const char &rhs) const;
  inline bool operator!=(const char &rhs) const;
  inline bool operator<(const char &rhs) const;
  friend bool operator==(const char &lhs, const MutationBase &rhs);
  friend bool operator!=(const char &lhs, const MutationBase &rhs);
  friend bool operator<(const char &lhs, const MutationBase &rhs);

  friend MutationBase::BitArray operator&(const BitArray &lhs, const BitArray &rhs);

  struct DNA {
    static constexpr size_t DNACount = 4;
    static const MutationBase A, C, G, T, N;
    static const char ambiguous_char;
    static const std::array<MutationBase, DNACount> bases;
    static const std::map<MutationBase, char> mut_to_char_map;
    static const std::map<MutationBase, MutationBase> complement_map;
  };

  BitArray value = {0, 0, 0, 0};
};

namespace std {
template <>
struct hash<MutationBase> {
  std::size_t operator()(const MutationBase &obj) const { return obj.ToChar(); }
};
}  // namespace std

#include "larch/impl/madag/mutation_base_impl.hpp"
