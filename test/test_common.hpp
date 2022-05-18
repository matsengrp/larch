#pragma once

#include <functional>
#include <string>
#include <stdexcept>

struct Test {
  std::function<void()> entry;
  std::string name;
};

bool add_test(const Test& test) noexcept;

inline void assert_true(bool expr, const std::string& what) {
  if (!expr) throw std::runtime_error(what);
}

inline void assert_false(bool expr, const std::string& what) {
  if (expr) throw std::runtime_error(what);
}

template <typename L, typename R>
inline void assert_equal(L&& l, R&& r, const std::string& what) {
  if (l != r) throw std::runtime_error(what);
}
