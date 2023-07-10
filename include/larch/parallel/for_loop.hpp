#pragma once

#include "larch/parallel/scheduler.hpp"

Scheduler& DefaultScheduler();

template <typename Lambda>
void seq_for_each(size_t size, Lambda&& lambda) {
  for (size_t i = 0; i < size; ++i) {
    lambda(i, 0);
  }
}

template <typename Lambda>
void parallel_for_each(size_t size, Lambda&& lambda) {
#if 0
  std::vector<std::thread> workers;
  std::atomic<size_t> iteration{0};
  for (size_t i = 0; i < std::thread::hardware_concurrency(); ++i) {
    workers.push_back(std::thread([&, i] {
      while (true) {
        size_t iter = iteration.fetch_add(1);
        if (iter >= size) {
          break;
        }
        lambda(iter, i);
      }
    }));
  }
  for (auto& i : workers) {
    i.join();
  }
#else
  Task task(DefaultScheduler(), [&](size_t i, size_t worker) {
    if (i >= size) {
      return false;
    }
    lambda(i, worker);
    return true;
  });
  DefaultScheduler().AddTask(task);
  task.Join(DefaultScheduler());
#endif
}
