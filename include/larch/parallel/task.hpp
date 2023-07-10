#pragma once

#include <atomic>
#include <mutex>

class Scheduler;

class TaskBase {
 public:
  TaskBase(const TaskBase&) = delete;
  TaskBase(TaskBase&&) = delete;
  TaskBase& operator=(const TaskBase&) = delete;
  TaskBase& operator=(TaskBase&&) = delete;
  inline explicit TaskBase(Scheduler& scheduler);
  virtual ~TaskBase() {}

  inline size_t GetId() const;

 protected:
  friend class Scheduler;
  virtual bool Run(size_t worker) = 0;
  virtual bool Finish(size_t workers_count) = 0;

 private:
  const size_t id_;
};

inline bool operator<(const std::reference_wrapper<TaskBase>& lhs,
                      const std::reference_wrapper<TaskBase>& rhs) {
  return std::addressof(lhs.get()) < std::addressof(rhs.get());
}

template <typename F>
class Task : public TaskBase {
 public:
  Task(Scheduler& scheduler, F&& func);
  ~Task();

  void Join();                      // Join by waiting
  void Join(Scheduler& scheduler);  // Join by working

 private:
  bool Run(size_t worker) override;
  bool Finish(size_t workers_count) override;

  F func_;
  std::atomic<size_t> iteration_ = 0;
  std::atomic<size_t> finished_ = 0;
  bool done_ = false;
  std::recursive_mutex done_mtx_;
  std::condition_variable_any is_done_;
};

