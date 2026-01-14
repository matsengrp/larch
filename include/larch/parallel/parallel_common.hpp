#pragma once

#include "larch/common.hpp"

#include <execution>
#include <latch>

#ifndef DISABLE_PARALLELISM
#include <taskflow/taskflow.hpp>
#endif

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
void SeqForEach(Range&& range, F&& func) {
  ranges::for_each(std::forward<Range>(range), std::forward<F>(func));
}

#ifndef DISABLE_PARALLELISM
inline tf::Executor& GetTaskflowExecutor() {
  static tf::Executor executor;
  return executor;
}
#endif

template <typename Range, typename F>
void ParallelForEach(Range&& range, F&& func) {
#ifdef DISABLE_PARALLELISM
  SeqForEach(std::forward<Range>(range), std::forward<F>(func));
#else
  std::vector vec = ranges::to_vector(range);
  if (vec.empty()) {
    return;
  }
  std::latch done{static_cast<std::ptrdiff_t>(vec.size())};
  auto& executor = GetTaskflowExecutor();
  for (auto& item : vec) {
    executor.silent_async([&func, &item, &done]() {
      func(item);
      done.count_down();
    });
  }
  done.wait();
#endif
}
