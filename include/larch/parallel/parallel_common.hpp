#pragma once

#include <mutex>
#include <shared_mutex>

#include <tbb/parallel_for_each.h>

#include "larch/common.hpp"

template <typename M>
auto ReadLock(M& mutex) {
  if constexpr (std::is_same_v<M, std::mutex>) {
    return std::unique_lock<M>{mutex};
  } else {
    return std::shared_lock<M>{mutex};
  }
}

template <typename M>
std::unique_lock<M> WriteLock(M& mutex) {
  return std::unique_lock{mutex};
}

template <typename T, typename M = std::shared_mutex>
class SharedState {
 public:
  SharedState() = default;
  explicit SharedState(T&& target) : target_{std::forward<T>(target)} {}

  template <typename F, typename... Args>
  decltype(auto) Read(F&& func, Args&&... args) const {
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
  decltype(auto) Write(F&& func, Args&&... args) {
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

template <typename Range, typename F>
void ParallelForEach(Range&& range, F&& func) {
  // tbb::parallel_for_each(std::forward<Range>(range), std::forward<F>(func));
  // ranges::for_each(std::forward<Range>(range), std::forward<F>(func));
  std::atomic<size_t> iter{0};
  const size_t count = size_t{range.size()};
  std::vector<std::thread> workers;
  for (size_t i = 0; i < 32; ++i) {
    workers.emplace_back(std::thread([&] {
      while (true) {
        const size_t it = iter.fetch_add(1);
        if (it >= count) {
          break;
        }
        func(range.at(it));
      }
    }));
  }
  for (auto& i : workers) {
    i.join();
  }
}
