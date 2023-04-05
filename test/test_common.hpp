#pragma once

#include <functional>
#include <string>
#include <stdexcept>
#include <iostream>

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
  if (not(l == r)) throw std::runtime_error(what);
}

inline bool test_true(bool expr, const std::string& what) {
  if (!expr) {
    std::cout << "TEST_FAILURE: " << what << std::endl;
  }
  return expr;
}

inline bool test_false(bool expr, const std::string& what) {
  if (expr) {
    std::cout << "TEST_FAILURE: " << what << std::endl;
  }
  return !expr;
}
