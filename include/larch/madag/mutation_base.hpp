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
  using BitArray = std::array<bool, 2>;

  inline MutationBase(const BitArray m_value);
  inline MutationBase(const char m_char_in);

  inline MutationBase GetComplementaryBase() const;

  inline char ToChar() const;
  static inline std::string ToString(std::vector<MutationBase> m_in);
  friend std::ostream &operator<<(std::ostream &os, const MutationBase m_in);
  friend std::ostream &operator<<(std::ostream &os,
                                  const std::vector<MutationBase> m_in);

  inline bool operator==(const MutationBase rhs) const;
  inline bool operator<(const MutationBase rhs) const;
  inline bool operator==(const BitArray rhs) const;
  inline bool operator==(const char rhs) const;
  inline bool operator!=(const char rhs) const;

  struct DNA {
    static constexpr size_t DNACount = 4;
    static const MutationBase A, C, G, T;
    static const std::array<MutationBase, DNACount> bases;
    static const std::map<MutationBase, char> mut_to_char_map;
    static const std::map<MutationBase, MutationBase> complement_map;
  };

  BitArray value;
};

#include "larch/impl/madag/mutation_base_impl.hpp"
