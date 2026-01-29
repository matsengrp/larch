#pragma once

#include <netam/include-fmt.hpp>
#include "larch/common.hpp"

#include <iostream>
#include <source_location>
#include <limits>

namespace netam {

[[noreturn]] inline void fail(std::string_view fmt_str, auto&&... args) {
  fmt::println("Failed: {}", fmt::vformat(fmt_str, fmt::make_format_args(args...)));
  std::cout << std::flush;
  throw std::runtime_error{fmt::vformat(fmt_str, fmt::make_format_args(args...))};
}

template <typename To, typename From>
  requires(std::integral<From> and std::integral<To>)
constexpr To checked_cast(From x) {
  if constexpr (std::is_signed_v<From> == std::is_signed_v<To>) {
    if (x < std::numeric_limits<To>::min()) {
      goto err;
    }
    if (x > std::numeric_limits<To>::max()) {
      goto err;
    }
  } else if constexpr (std::is_unsigned_v<To>) {
    if (x < 0) {
      goto err;
    }
    // TODO
  } else {
    // TODO
  }
  return static_cast<To>(x);
err:
  fail("cast out of bounds");
}

template <typename To, typename From>
  requires((sizeof(To) < sizeof(From)) and
           (std::is_signed_v<From> == std::is_signed_v<To>))
constexpr To narrowing_cast(From x) {
  return checked_cast<To>(x);
}

template <typename From>
  requires(std::is_unsigned_v<From>)
constexpr std::make_signed_t<From> signed_cast(From x) {
  return checked_cast<std::make_signed_t<From>>(x);
}

template <typename From>
  requires(std::is_signed_v<From>)
constexpr std::make_unsigned_t<From> unsigned_cast(From x) {
  return checked_cast<std::make_unsigned_t<From>>(x);
}

}  // namespace netam
