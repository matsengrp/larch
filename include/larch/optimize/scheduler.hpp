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
  virtual void Run(size_t worker) = 0;
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
  void Run(size_t worker) override;
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
  inline explicit Scheduler(size_t workers_count = std::thread::hardware_concurrency());
  inline ~Scheduler();

  Scheduler(const Scheduler&) = delete;
  Scheduler(Scheduler&&) = delete;
  Scheduler& operator=(const Scheduler&) = delete;
  Scheduler& operator=(Scheduler&&) = delete;

  inline void AddTask(TaskBase& task);

  inline size_t WorkersCount() const;

 private:
  static inline void Worker(Scheduler& self, size_t id);
  const size_t workers_count_;
  std::vector<std::thread> workers_;
  std::atomic<bool> destroy_ = false;
  std::mutex mtx_;
  std::deque<std::reference_wrapper<TaskBase>> queue_;
  std::set<std::reference_wrapper<TaskBase>> finished_tasks_;
  std::condition_variable queue_not_empty_;
};

template <typename T>
class Reduction {
 public:
  explicit Reduction(size_t workers_count) : size_{0} { data_.resize(workers_count); }

  template <typename... Args>
  T& Emplace(size_t worker, Args&&... args) {
    size_.fetch_add(1);
    return data_.at(worker).emplace_back(std::forward<Args>(args)...);
  }

  auto GetAll() { return data_ | ranges::views::join; }
  auto GetAll() const { return data_ | ranges::views::join; }

  size_t Size() const { return size_.load(); }

 private:
  std::vector<std::deque<T>> data_;
  std::atomic<size_t> size_;
};

#include "larch/impl/optimize/scheduler_impl.hpp"
