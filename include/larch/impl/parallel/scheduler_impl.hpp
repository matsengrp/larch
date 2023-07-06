
TaskBase::TaskBase(Scheduler& scheduler) : id_{scheduler.NewTaskId()} {}

size_t TaskBase::GetId() const { return id_; }

template <typename F>
Task<F>::Task(Scheduler& scheduler, F&& func)
    : TaskBase{scheduler}, func_{std::forward<F>(func)} {}

template <typename F>
Task<F>::~Task() {
  Join();
}

template <typename F>
void Task<F>::Join() {
  std::unique_lock lock{done_mtx_};
  while (not done_) {
    is_done_.wait(lock);
  }
}

template <typename F>
void Task<F>::Join(Scheduler& scheduler) {
  if (not scheduler.JoinTask(*this)) {
    Join();
  }
}

template <typename F>
bool Task<F>::Run(size_t worker) {
  return func_(iteration_.fetch_add(1), worker);
}

template <typename F>
bool Task<F>::Finish(size_t workers_count) {
  if (finished_.fetch_add(1) + 1 == workers_count) {
    std::unique_lock lock{done_mtx_};
    done_ = true;
    is_done_.notify_all();
    return true;
  }
  return false;
}

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
  for (size_t id = 0; id < workers_count_; ++id) {
    Worker& worker = workers_.at(id);
    std::unique_lock lock{worker.mtx};
    if (worker.thread == nullptr) {
      throw std::runtime_error("Fail");
    }
    if (worker.thread->get_id() == std::this_thread::get_id()) {
      lock.unlock();
      WorkUntil(
          id,
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
    }
  }
  return false;
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

template <typename T, typename WorkerId>
Reduction<T, WorkerId>::Reduction(size_t workers_count) : size_{0} {
  if constexpr (UseVector) {
    data_.resize(workers_count);
  }
}

template <typename T, typename WorkerId>
Reduction<T, WorkerId>::Reduction() : size_{0} {
  static_assert(not UseVector);
}

template <typename T, typename WorkerId>
template <typename... Args>
T& Reduction<T, WorkerId>::Emplace(WorkerId worker, Args&&... args) {
  size_.fetch_add(1);
  if constexpr (UseVector) {
    return data_[worker].emplace_back(std::forward<Args>(args)...);
  } else {
    return data_.At(worker).emplace_back(std::forward<Args>(args)...);
  }
}

template <typename T, typename WorkerId>
auto Reduction<T, WorkerId>::Get() {
  if constexpr (UseVector) {
    return data_ | ranges::views::all;
  } else {
    return data_.All() | ranges::views::values;
  }
}

template <typename T, typename WorkerId>
auto Reduction<T, WorkerId>::GetAll() {
  return SizedView{Get() | ranges::views::join, Size()};
}

template <typename T, typename WorkerId>
size_t Reduction<T, WorkerId>::Size() const {
  return size_.load();
}

template <typename T, typename WorkerId>
bool Reduction<T, WorkerId>::Empty() const {
  return size_.load() > 0;
}

template <typename T, typename WorkerId>
size_t Reduction<T, WorkerId>::WorkersCount() const {
  return data_.size();
}

template <typename T, typename WorkerId>
void Reduction<T, WorkerId>::Clear() {
  size_.store(0);
  data_.Clear();
}

template <typename T>
template <typename... Args>
Snapshot<T>::Snapshot(Args&&... args) : data_{new T{std::forward<Args>(args)...}} {}

template <typename T>
Snapshot<T>::~Snapshot() {
  delete data_.load();
}

template <typename T>
T& Snapshot<T>::Get() {
  return *data_.load();
}

template <typename T>
template <typename Lambda, typename... Args>
void Snapshot<T>::Take(Lambda&& lambda, Args&&... args) {
  T* next = new T{std::forward<Args>(args)...};
  std::unique_ptr<T> curr{std::atomic_exchange(&data_, next)};
  lambda(*curr);
}
