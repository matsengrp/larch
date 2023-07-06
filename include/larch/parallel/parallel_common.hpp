#pragma once

#include <mutex>

#include "larch/common.hpp"

template <typename T>
class Accessor {
 public:
  Accessor(T& value, std::mutex& mtx) : value_{value}, mtx_{mtx} {}

  template <typename Lambda>
  void Get(Lambda&& lambda) {
    std::unique_lock lock{mtx_};
    lambda(value_);
  }

  template <typename Lambda, typename... Args>
  decltype(auto) Get2(Lambda&& lambda, Args&&... args) {
    std::unique_lock lock{mtx_};
    return lambda(value_, std::forward<Args>(args)...);
  }

 private:
  T& value_;
  std::mutex& mtx_;
};
