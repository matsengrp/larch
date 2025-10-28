#pragma once

#include "larch/common.hpp"

#include <execution>

// #include <tbb/parallel_for_each.h>

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
  std::vector vec = ranges::to_vector(range);
  std::for_each(std::execution::par, 
    std::begin(vec), std::end(vec), std::forward<F>(func));
  // tbb::parallel_for_each(std::begin(vec), std::end(vec),
  //               std::forward<F>(func));
}

template <typename Range, typename F>
void SeqForEach(Range&& range, F&& func) {
  ranges::for_each(std::forward<Range>(range), std::forward<F>(func));
}
