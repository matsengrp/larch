Scheduler::Scheduler(size_t workers_count)
    : workers_count_{workers_count}, workers_{workers_count} {
  for (size_t i = 0; i < workers_count; ++i) {
    std::unique_lock lock{workers_.at(i).mtx};
    workers_.at(i).thread =
        std::make_unique<std::thread>(&Scheduler::WorkerThread, this, i);
  }
}

Scheduler::~Scheduler() {
  destroy_.store(true);
  for (Worker& worker : workers_) {
    std::unique_lock lock{worker.mtx};
    worker.has_work.notify_all();
  }
  for (auto& i : workers_) {
    i.thread->join();
  }
}

void Scheduler::AddTask(TaskBase& task) {
  {
    std::unique_lock lock{running_tasks_mtx_};
    running_tasks_.insert(task.GetId());
  }
  for (Worker& worker : workers_) {
    std::unique_lock lock{worker.mtx};
    worker.queue.push_back({std::ref(task), task.GetId()});
    worker.has_work.notify_all();
  }
}

size_t Scheduler::WorkersCount() const { return workers_count_; }

bool Scheduler::JoinTask(TaskBase& task) {
  auto worker_id = FindWorkerId();
  if (worker_id.has_value()) {
    WorkUntil(
        worker_id.value(),
        [&] {
          if (destroy_.load()) {
            return true;
          }
          std::unique_lock tasks_lock{running_tasks_mtx_};
          return running_tasks_.find(task.GetId()) == running_tasks_.end();
        },
        [&](QueueItem& item) {
          return not item.done and item.task.get().GetId() != task.GetId();
        });
    return true;
  } else {
    return false;
  }
}

size_t Scheduler::NewTaskId() { return task_ids_.fetch_add(1); }

void Scheduler::WorkerThread(const size_t id) {
  WorkUntil(
      id, [&] { return destroy_.load(); },
      [&](QueueItem& item) { return not item.done; });
}

template <typename Until, typename Accept>
void Scheduler::WorkUntil(size_t id, Until&& until, Accept&& accept) {
  Worker& worker = workers_.at(id);
  std::unique_lock lock{worker.mtx};
  while (not until()) {
    while (worker.queue.empty()) {
      worker.has_work.wait(lock);
      if (until()) {
        return;
      }
    }
    while (not worker.queue.empty() and worker.queue.front().done) {
      worker.queue.pop_front();
    }
    if (worker.queue.empty()) {
      continue;
    }
    for (QueueItem& item : worker.queue) {
      if (not accept(item)) {
        continue;
      }
      const size_t task_id = item.task.get().GetId();
      lock.unlock();
      if (not item.task.get().Run(id)) {
        {
          std::unique_lock tasks_lock{running_tasks_mtx_};
          if (item.task.get().Finish(workers_count_)) {
            running_tasks_.erase(task_id);
          }
        }
        lock.lock();
        item.done = true;
      } else {
        lock.lock();
      }
      break;
    }
  }
}

std::optional<size_t> Scheduler::FindWorkerId() const {
  for (size_t id = 0; id < workers_count_; ++id) {
    const Worker& worker = workers_.at(id);
    std::unique_lock lock{worker.mtx};
    if (worker.thread == nullptr) {
      throw std::runtime_error("Fail");
    }
    if (worker.thread->get_id() == std::this_thread::get_id()) {
      return id;
    }
  }
  return std::nullopt;
}
