#pragma once

#include <mutex>
#include <shared_mutex>

#include "larch/common.hpp"

template <typename M>
using SharedLock = std::conditional_t<std::is_same_v<M, std::shared_mutex>,
                                      std::shared_lock<M>, std::unique_lock<M>>;

template <typename T, typename M>
class Accessor {
 public:
  Accessor(T& value, M& mtx) : value_{value}, mtx_{mtx} {}

  template <typename Lambda, typename... Args>
  decltype(auto) GetExclusive(Lambda&& lambda, Args&&... args) {
    std::unique_lock lock{mtx_};
    if constexpr (std::is_void_v<decltype(lambda(value_,
                                                 std::forward<Args>(args)...))>) {
      lambda(value_, std::forward<Args>(args)...);
    } else {
      return lambda(value_, std::forward<Args>(args)...);
    }
  }

  template <typename Lambda, typename... Args>
  decltype(auto) GetShared(Lambda&& lambda, Args&&... args) {
    SharedLock<M> lock{mtx_};
    if constexpr (std::is_void_v<decltype(lambda(value_,
                                                 std::forward<Args>(args)...))>) {
      lambda(value_, std::forward<Args>(args)...);
    } else {
      return lambda(value_, std::forward<Args>(args)...);
    }
  }

 private:
  T& value_;
  M& mtx_;
};
