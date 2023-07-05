
template <typename F>
Task<F>::Task(F&& func) : func_{std::forward<F>(func)} {}

template <typename F>
void Task<F>::Join(Scheduler& scheduler) {
  if (not scheduler.WorkUntilDone(*this)) {
    while (not done_.load()) {
      std::unique_lock lock{done_mtx_};
      is_done_.wait(lock);
    }
  }
}

template <typename F>
bool Task<F>::Run(size_t worker) {
  return func_(iteration_.fetch_add(1), worker);
}

template <typename F>
void Task<F>::Finish(size_t workers_count) {
  if (finished_.fetch_add(1) + 1 == workers_count) {
    done_.store(true);
    std::unique_lock lock{done_mtx_};
    is_done_.notify_all();
  }
}

template <typename F>
bool Task<F>::IsDone() const {
  return done_.load();
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
  for (auto& i : workers_) {
    i.thread->join();
  }
}

size_t Scheduler::AddTask(TaskBase& task) {
  size_t task_id = task_ids_.fetch_add(1);
  for (Worker& worker : workers_) {
    std::unique_lock lock{worker.mtx};
    worker.queue.push_back({std::ref(task), task_id});
    worker.has_work.notify_all();
  }
  return task_id;
}

size_t Scheduler::WorkersCount() const { return workers_count_; }

bool Scheduler::WorkUntilDone(TaskBase& task) {
  for (size_t id = 0; id < workers_count_; ++id) {
    Worker& worker = workers_.at(id);
    std::unique_lock lock{worker.mtx};
    if (worker.thread == nullptr) {
      throw std::runtime_error("Fail");
    }
    if (worker.thread->get_id() == std::this_thread::get_id()) {
      lock.unlock();
      WorkUntil(id, [&] { return task.IsDone(); });
      return true;
    }
  }
  return false;
}

void Scheduler::WorkerThread(const size_t id) {
  WorkUntil(id, [&] { return destroy_.load(); });
}

template <typename Lambda>
void Scheduler::WorkUntil(size_t id, Lambda&& until) {
  Worker& worker = workers_.at(id);
  std::unique_lock lock{worker.mtx};
  while (not until()) {
    while (worker.queue.empty()) {
      worker.has_work.wait(lock);
    }
    while (not worker.queue.empty() and worker.queue.front().done) {
      worker.queue.pop_front();
    }
    if (worker.queue.empty()) {
      continue;
    }
    for (QueueItem& item : worker.queue) {
      if (item.done) {
        continue;
      }
      lock.unlock();
      if (not item.task.get().Run(id)) {
        item.task.get().Finish(workers_count_);
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
#ifdef USE_TSAN
  std::unique_lock lock{tsan_mtx_};
#endif
  return data_[worker].emplace_back(std::forward<Args>(args)...);
}

template <typename T, typename WorkerId>
auto Reduction<T, WorkerId>::Get() {
  if constexpr (UseVector) {
    return data_ | ranges::views::all;
  } else {
    return data_ | ranges::views::values;
  }
}

template <typename T, typename WorkerId>
auto Reduction<T, WorkerId>::GetAll() {
#ifdef USE_TSAN
  std::unique_lock lock{tsan_mtx_};
#endif
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
#ifdef USE_TSAN
  std::unique_lock lock{tsan_mtx_};
#endif
  data_.clear();
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
