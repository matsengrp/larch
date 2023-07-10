#pragma once

#include <thread>
#include <condition_variable>
#include <deque>
#include <memory>
#include <functional>
#include <set>
#include <type_traits>
#include <optional>

#include "larch/parallel/node_hashset.hpp"
#include "larch/parallel/task.hpp"

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
    mutable std::mutex mtx = {};
    std::condition_variable has_work = {};
  };

  inline size_t NewTaskId();

  inline void WorkerThread(size_t id);

  template <typename Until, typename Accept>
  void WorkUntil(size_t id, Until&& until, Accept&& accept);

  inline std::optional<size_t> FindWorkerId() const;

  const size_t workers_count_;
  std::vector<Worker> workers_;
  std::atomic<bool> destroy_ = false;
  std::atomic<size_t> task_ids_ = 0;

  std::set<size_t> running_tasks_;
  std::mutex running_tasks_mtx_;  // TODO concurrent flat set
};

#include "larch/impl/parallel/task_impl.hpp"
#include "larch/impl/parallel/scheduler_impl.hpp"
