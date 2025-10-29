#pragma once

#include "larch/common.hpp"

#include <execution>

#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>

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
  // std::for_each(std::execution::par,
  //   std::begin(vec), std::end(vec), std::forward<F>(func));
  using block = tbb::blocked_range<size_t>;
  tbb::parallel_for(block{0, vec.size()}, [&vec, func](const block& x) {
    for (size_t i = x.begin(); i != x.end(); ++i) {
      func(vec.at(i));
    }
  });
}

template <typename Range, typename F>
void SeqForEach(Range&& range, F&& func) {
  ranges::for_each(std::forward<Range>(range), std::forward<F>(func));
}
