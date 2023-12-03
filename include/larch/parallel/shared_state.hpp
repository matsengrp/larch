#pragma once

#include "larch/parallel/rw_lock.hpp"

template <typename T, typename M = std::shared_mutex>
class SharedState {
 public:
  SharedState() = default;

  SharedState(SharedState&& other)
      : target_{[](SharedState&& o) {
          auto lock = ReadLock(o.mutex_);
          return std::move(o.target_);
        }(std::move(other))} {}
  explicit SharedState(T&& target) : target_{std::forward<T>(target)} {}

  SharedState& operator=(SharedState&& other) {
    auto wl = WriteLock(mutex_);
    auto rl = ReadLock(other.mutex_);
    target_ = std::move(other.target_);
    return *this;
  };

  template <typename F, typename... Args>
  decltype(auto) ReadShared(F&& func, Args&&... args) const {
    auto lock = ReadLock(mutex_);
    if constexpr (std::is_void_v<decltype(std::invoke(std::forward<F>(func),
                                                      static_cast<const T&>(target_),
                                                      std::forward<Args>(args)...))>) {
      std::invoke(std::forward<F>(func), static_cast<const T&>(target_),
                  std::forward<Args>(args)...);
    } else {
      return std::invoke(std::forward<F>(func), static_cast<const T&>(target_),
                         std::forward<Args>(args)...);
    }
  }

  template <typename F, typename... Args>
  decltype(auto) WriteShared(F&& func, Args&&... args) {
    auto lock = WriteLock(mutex_);
    if constexpr (std::is_void_v<decltype(std::invoke(std::forward<F>(func),
                                                      static_cast<T&>(target_),
                                                      std::forward<Args>(args)...))>) {
      std::invoke(std::forward<F>(func), static_cast<T&>(target_),
                  std::forward<Args>(args)...);
    } else {
      return std::invoke(std::forward<F>(func), static_cast<T&>(target_),
                         std::forward<Args>(args)...);
    }
  }

 private:
  mutable M mutex_;
  T target_;
};
