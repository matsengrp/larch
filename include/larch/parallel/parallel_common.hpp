#pragma once

#include "larch/common.hpp"

#include <execution>
#include <latch>

#ifndef DISABLE_PARALLELISM
#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>
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
  auto& executor = GetTaskflowExecutor();
  if (executor.this_worker_id() >= 0) {
    tf::Taskflow taskflow;
    taskflow.for_each(vec.begin(), vec.end(), [&func](auto& item) { func(item); });
    executor.corun(taskflow);
  } else {
    std::atomic<size_t> counter{vec.size()};
    for (auto& item : vec) {
      executor.silent_dependent_async([&func, &item, &counter]() {
        func(item);
        counter.fetch_sub(1, std::memory_order_release);
      });
    }
    while (counter.load(std::memory_order_acquire) > 0) {
      std::this_thread::yield();
    }
  }
#endif
}
