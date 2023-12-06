#pragma once

#include <mutex>
#include <shared_mutex>
#include <execution>
#include <thread>
#include <atomic>

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

template <typename Range, typename F>
void ParallelForEach(Range&& range, F&& func) {
  // ranges::for_each(std::forward<Range>(range), std::forward<F>(func));
  std::for_each(std::execution::par_unseq, std::begin(range), std::end(range),
                std::forward<F>(func));
}
