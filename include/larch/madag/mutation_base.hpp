/**
 * MutationBase represents a DNA base value in the MADAG. Does not support ambiguity or
 * gap characters. Bases are encoded as a two-bit array. Bases are ordered such that
 * complementary bases are represented by inverting the bits.
 */

#pragma once

#include <string_view>
#include <vector>

#include "larch/contiguous_map.hpp"

struct MutationBase {
  using type = std::uint32_t;

  MutationBase() : value{0} {}  // TODO

  constexpr MutationBase(char x) noexcept : value{x == 'N' ? ~zero : mask(x)} {}

  template <typename... Chars>
  constexpr MutationBase(Chars... x) noexcept : value{zero | (mask(x) | ...)} {
    static_assert(sizeof...(Chars) > 1);
    static_assert((std::is_same_v<Chars, char> and ...));
  }

  constexpr char ToChar() const noexcept {
    const int ctz = __builtin_ctz(value);
    const type bit = one << ctz;
    if ((value & bit) == value) {
      return 'A' + static_cast<char>(ctz);
    } else {
      return 'N';
    }
  }

  inline bool IsCompatible(const MutationBase& rhs) const;
  inline bool IsAmbiguous() const;

  inline MutationBase GetComplementaryBase() const;
  inline MutationBase GetFirstBase() const;
  inline MutationBase GetCommonBases(const MutationBase& rhs) const;
  inline MutationBase GetFirstCommonBase(const MutationBase& rhs) const;

  static inline std::string ToString(const std::vector<MutationBase>& m_in);

  inline bool operator==(const MutationBase& rhs) const;
  inline bool operator!=(const MutationBase& rhs) const;
  inline bool operator<(const MutationBase& rhs) const;
  inline bool operator==(const char& rhs) const;
  inline bool operator!=(const char& rhs) const;
  inline bool operator<(const char& rhs) const;
  friend bool operator==(const char& lhs, const MutationBase& rhs);
  friend bool operator!=(const char& lhs, const MutationBase& rhs);
  friend bool operator<(const char& lhs, const MutationBase& rhs);

 private:
  MutationBase(type x) : value{x} {}

  constexpr int Count() const noexcept { return __builtin_popcount(value); }

  constexpr bool Test(char x) const noexcept { return (value & mask(x)) != 0; }

  constexpr void Set(char x) noexcept { value |= mask(x); }

  constexpr void Clear(char x) noexcept { value &= ~mask(x); }

  constexpr void SetOnly(char x) noexcept { value = mask(x); }

  constexpr void SetAll() noexcept { value = ~zero; }

  constexpr std::optional<MutationBase> Common(MutationBase rhs) const noexcept {
    type result = value & rhs.value;
    if (result == 0) {
      return std::nullopt;
    }
    return MutationBase{result};
  }

  constexpr char GetFirst() const noexcept {
    return 'A' + static_cast<char>(__builtin_ctz(value));
  }

  static constexpr type mask(char x) noexcept { return one << (x - 'A'); }

  type value;
  static constexpr const type zero = 0;
  static constexpr const type one = 1;
};

static_assert(std::is_trivially_copyable_v<MutationBase>);

namespace std {
template <>
struct hash<MutationBase> {
  std::size_t operator()(const MutationBase& obj) const {
    return static_cast<std::size_t>(obj.ToChar());
  }
};
}  // namespace std

#include "larch/impl/madag/mutation_base_impl.hpp"
