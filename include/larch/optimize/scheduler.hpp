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
#include <map>
#include <type_traits>

class TaskBase {
 public:
  virtual ~TaskBase() {}

 protected:
  friend class Scheduler;
  virtual bool Run(size_t worker) = 0;
  virtual void Finish() = 0;
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
  bool Run(size_t worker) override;
  void Finish() override;

  F func_;
  std::atomic<size_t> iteration_;
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
  using QueueItem = std::pair<std::reference_wrapper<TaskBase>, size_t>;

  inline void Worker(size_t id);
  inline bool IsTaskRunning(size_t task_id) const;

  const size_t workers_count_;
  std::vector<std::thread> workers_;
  std::vector<size_t> running_tasks_;
  std::atomic<bool> destroy_ = false;
  std::atomic<size_t> ids_ = 0;
  std::mutex mtx_;
  std::deque<QueueItem> queue_;
  std::condition_variable queue_not_empty_;
  std::set<size_t> done_tasks_;
};

template <typename V>
class SizedView {
 public:
  SizedView(V&& view, size_t size) : view_{view}, size_{size} {}
  decltype(auto) begin() { return view_.begin(); }
  decltype(auto) end() { return view_.end(); }
  size_t size() const { return size_; }

 private:
  V view_;
  size_t size_;
};

template <typename T, typename WorkerId = size_t>
class Reduction {
  static constexpr bool UseVector = std::is_same_v<WorkerId, size_t>;
  using Container = std::conditional_t<UseVector, std::vector<std::deque<T>>,
                                       ConcurrentUnorderedMap<WorkerId, std::deque<T>>>;

 public:
  explicit Reduction(size_t workers_count);

  Reduction();

  template <typename... Args>
  T& Emplace(WorkerId worker, Args&&... args);

  auto Get();

  auto GetAll();

  size_t Size() const;
  bool Empty() const;
  size_t WorkersCount() const;

  void Clear();

 private:
  Container data_;
  std::atomic<size_t> size_;
#ifdef USE_TSAN
  std::mutex mtx_;
#endif
};

template <typename T>
class Snapshot {
 public:
  template <typename... Args>
  explicit Snapshot(Args&&... args);

  ~Snapshot();

  T& Get();

  template <typename Lambda, typename... Args>
  void Take(Lambda&& lambda, Args&&... args);

 private:
  std::atomic<T*> data_;
};

#include "larch/impl/optimize/scheduler_impl.hpp"
