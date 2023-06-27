#pragma once

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <memory>
#include <atomic>
#include <functional>
#include <iostream>

class TaskBase {
 public:
  virtual ~TaskBase() {}

 protected:
  friend class Scheduler;
  virtual void Run() = 0;
  virtual bool CanIterate() const = 0;
  virtual void Finish() = 0;
};

template <typename F>
class Task : public TaskBase {
 public:
  Task(F&& func) : func_{std::forward<F>(func)} {}

  void Join() {
    std::unique_lock lock{mtx_};
    while (not is_finished_) {
      finished_.wait(lock);
    }
  }

 private:
  void Run() override {
    if (not func_(iteration_.fetch_add(1))) {
      can_iterate_.store(false);
    }
  }

  bool CanIterate() const override { return can_iterate_.load(); }

  void Finish() override {
    std::unique_lock lock{mtx_};
    if (not is_finished_) {
      is_finished_ = true;
      finished_.notify_all();
    }
  }

  F func_;
  std::atomic<size_t> iteration_ = 0;
  std::atomic<bool> can_iterate_ = true;
  std::mutex mtx_;
  std::condition_variable finished_;
  bool is_finished_ = false;
};

class Scheduler {
 public:
  inline ~Scheduler();

  inline void Start();

  inline void AddTask(TaskBase& task);

 private:
  static inline void Worker(Scheduler& self);
  std::vector<std::thread> workers_;
  std::atomic<bool> destroy_ = false;
  std::mutex mtx_;
  std::deque<std::reference_wrapper<TaskBase>> queue_;
  std::condition_variable queue_not_empty_;
};

#include "larch/impl/optimize/scheduler_impl.hpp"
