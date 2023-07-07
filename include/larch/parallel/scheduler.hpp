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
#include <type_traits>
#include <optional>
#include <shared_mutex>

#include "larch/parallel/node_hashset.hpp"
#include "larch/parallel/node_hashmap.hpp"

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
  std::mutex done_mtx_;
  std::condition_variable is_done_;
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

  inline bool JoinTask(TaskBase& task);

 private:
  friend class TaskBase;

  struct QueueItem {
    std::reference_wrapper<TaskBase> task;
    const size_t id;
    bool done = false;
  };

  struct Worker {
    std::unique_ptr<std::thread> thread = nullptr;
    std::deque<QueueItem> queue = {};
    std::mutex mtx = {};
    std::condition_variable has_work = {};
  };

  inline size_t NewTaskId();

  inline void WorkerThread(size_t id);

  template <typename Until, typename Accept>
  void WorkUntil(size_t id, Until&& until, Accept&& accept);

  const size_t workers_count_;
  std::vector<Worker> workers_;
  std::atomic<bool> destroy_ = false;
  std::atomic<size_t> task_ids_ = 0;

  std::set<size_t> running_tasks_;
  std::mutex running_tasks_mtx_;  // TODO concurrent flat set
};

template <typename V>
class SizedView {
 public:
  SizedView(V&& view, size_t size) : view_{view}, size_{size} {}
  decltype(auto) begin() { return view_.begin(); }
  decltype(auto) end() { return view_.end(); }
  size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }

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
  ~Reduction();

  template <typename... Args>
  std::pair<T&, size_t> Emplace(WorkerId worker, Args&&... args);

  template <typename Lambda>
  void Consume(Lambda&& lambda);

  template <typename Lambda>
  void ConsumeBatches(Lambda&& lambda);

  size_t SizeApprox() const;
  size_t WorkersCount() const;

 private:
  struct Data {
    Container container;
    std::atomic<size_t> size = 0;
  };
  static auto GetRange(Data* data);
  static auto* MakeData();
  std::atomic<Data*> data_ = nullptr;
  std::shared_mutex mtx_;
};

#include "larch/impl/parallel/scheduler_impl.hpp"
