#pragma once

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <memory>
#include <atomic>
#include <functional>
#include <set>

class TaskBase {
 public:
  virtual ~TaskBase() {}

 protected:
  friend class Scheduler;
  virtual void Run() = 0;
  virtual bool CanIterate() const = 0;
  virtual bool Finish() = 0;
};

inline bool operator<(const std::reference_wrapper<TaskBase>& lhs,
                      const std::reference_wrapper<TaskBase>& rhs) {
  return std::addressof(lhs.get()) < std::addressof(rhs.get());
}

template <typename F>
class Task : public TaskBase {
 public:
  Task(F&& func);

  void Join();

 private:
  void Run() override;

  bool CanIterate() const override;

  bool Finish() override;

  F func_;
  std::atomic<size_t> iteration_;
  std::atomic<bool> can_iterate_;
  std::mutex mtx_;
  std::condition_variable finished_;
  bool is_finished_ = false;
};

class Scheduler {
 public:
  inline Scheduler();
  inline ~Scheduler();

  Scheduler(const Scheduler&) = delete;
  Scheduler(Scheduler&&) = delete;
  Scheduler& operator=(const Scheduler&) = delete;
  Scheduler& operator=(Scheduler&&) = delete;

  inline void AddTask(TaskBase& task);

 private:
  static inline void Worker(Scheduler& self);
  std::vector<std::thread> workers_;
  std::atomic<bool> destroy_ = false;
  std::mutex mtx_;
  std::deque<std::reference_wrapper<TaskBase>> queue_;
  std::set<std::reference_wrapper<TaskBase>> finished_tasks_;
  std::condition_variable queue_not_empty_;
};

#include "larch/impl/optimize/scheduler_impl.hpp"
